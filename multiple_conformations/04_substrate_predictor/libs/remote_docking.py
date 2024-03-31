#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging
import json
import time
import random
import shutil
import tarfile
import multiprocessing as mp

from libs import paths
from external import remote_processing as remote
from libs.docking import test_vina_results
from libs import utils
from libs.configuration import REMOTE, ENV
from libs.remote_docking_prepare import create_docking_batch_script

logger = logging.getLogger(__name__)


def candidates_remote_docking():
    vs = VirtualScreen()
    vs.set_up_screening()
    vs.start_screening()

    while not vs.check_jobs():
        time.sleep(REMOTE["sleep_time"])

    vs.sync_data_locally()


class VirtualScreen:
    def __init__(self):
        random.seed()
        self.job_list = []
        self.queue_manager = remote.QueueManager().get_handler()
        self.remote_root_dir = REMOTE["host"]["submit_dir"]
        self.dataname_local = REMOTE["docking_dataname"]
        self.datadir_local = os.path.join(paths.temp(), self.dataname_local)
        self.dataname_remote = self.dataname_local + "_" + str(random.randrange(int(REMOTE["batch_size"])))
        self.remote_screening_dir = os.path.join(self.remote_root_dir, self.dataname_remote)
        self.output_filename_pattern = self.queue_manager.get_output_filename_pattern()
        self.remote_job_name = "batch"
        self.source = os.path.join(self.datadir_local, self.dataname_local + ".tar.gz")
        self.dest = os.path.join(self.remote_screening_dir, self.dataname_local + ".tar.gz")
        self.remote_connection = remote.RemoteHandler()
        self._screening_status_storage = os.path.join(paths.temp(), "VS_storage.json")
        self.screening_finished = False
        self.failed_files = {}
        self.num_failed_jobs = 0
        self.current_jobs_statuses = {
            "F": 0,
            "R": 0,
            "Q": 0
        }
        self.names_to_delete = ["out.tar.gz", "log.tar.gz", self.output_filename_pattern]

    def __is_frequent_job_error(self):
        return True if self.num_failed_jobs >= len(self.job_list) / 2 else False

    def __remote_datadir_exists(self):
        return self.dataname_remote in self.remote_connection.listdir(self.remote_root_dir)

    def __delete_remote_screening_data(self):
        if self.screening_finished and self.__remote_datadir_exists():
            commands = []
            for job in self.job_list:
                if not (job.name in self.failed_files.keys() or self.remote_connection.file_info(job.location) is None):
                    logger.debug("Adding %s among folders for cleaning", job.location)
                    commands.append("rm -rf {}".format(job.location))

            if commands:
                self.remote_connection.execute_remote_command(commands)

            if self.__count_failed_files() == 0:
                self.remote_connection.remove_dir(self.remote_screening_dir)

    def __save_screening(self):
        tmp = {
            "jobs": self.__browse_jobs(),
            "dataname": self.dataname_remote,
            "num_failed_jobs": self.num_failed_jobs,
            "screening_finished": self.screening_finished
        }

        with open(self._screening_status_storage, "w") as out_stream:
            json.dump(tmp, out_stream)

        logger.debug("Screening saved to file %s", self._screening_status_storage)

    def __load_screening(self):
        with open(self._screening_status_storage, "r") as in_stream:
            jobs_data = json.load(in_stream)
            self.dataname_remote = jobs_data["dataname"]
            self.remote_screening_dir = os.path.join(self.remote_root_dir, self.dataname_remote)
            self.screening_finished = jobs_data["screening_finished"]
            self.num_failed_jobs = jobs_data["num_failed_jobs"]
            for job_data in jobs_data["jobs"]:
                job = RemoteJob(job_data["location"],  self.remote_connection, self.queue_manager, self.names_to_delete)
                if job_data["job_id"] is not None:
                    job.status = job_data["status"]
                    job.job_id = job_data["job_id"]
                self.job_list.append(job)

        if self.screening_finished:
            logger.info("Screening has already been performed.")
        elif self.__is_frequent_job_error():
            logger.error("More than half of the submitted jobs failed, terminating the screening")
        elif not self.__remote_datadir_exists():
            logger.warning("Screening directory %s of restarted process does not exist on remote host",
                           self.dataname_remote)
            logger.warning("Removing the restart file: %s and starting a new screening", self._screening_status_storage)
            utils.remove_existing(self._screening_status_storage)
            self.set_up_screening()

    def __parse_job_info(self, raw_input):
        job_info = {}
        if raw_input is None:
            return None
        elif raw_input == ['']:
            return ['']
        else:
            for line in raw_input:
                if str(self.remote_job_name) in line:
                    tmp = self.queue_manager.parse_job_info_line(line)
                    job_info[tmp["name"]] = tmp
            return job_info

    def __get_queue_status(self):
        logger.info("Asking queue manager on remote host for status of jobs.")
        commands = ["{}".format(self.queue_manager.get_queue_status_command())]
        raw_results = self.remote_connection.execute_remote_command(commands, True).decode('ascii').split("\n")
        status = self.__parse_job_info(raw_results)
        if len(self.job_list) != len(status):
            logger.debug("Queue manager status retrieved for %d/%d jobs only", len(status), len(self.job_list))
        return status

    def __evaluate_screening_status(self):
        if self.__is_frequent_job_error():
            logger.error("More than half of the submitted jobs failed, terminating the screening")
            self.__save_screening()
            utils.update_candidate_processing_status("docking-failed")
            return True

        if len(self.job_list) == self.current_jobs_statuses["F"]:
            logger.info("Screening finished")
            self.screening_finished = True
            self.__save_screening()
            return True
        else:
            logger.info("From %d jobs, %d are finished, %d are running, and %d are queued", len(self.job_list),
                        self.current_jobs_statuses["F"], self.current_jobs_statuses["R"],
                        self.current_jobs_statuses["Q"])
            self.__save_screening()
            return False

    def __browse_jobs(self):
        reports = []
        for job in self.job_list:
            reports.append(job.report())
        return reports

    def __resubmit_failed_job(self, job):
        logger.info("Verification of results failed for %s. \n Resubmitting this batch.", job.name)
        logger.info("This is %d failed job of all %d jobs.", self.num_failed_jobs, len(self.job_list))
        job.clean_job_dir()
        job.job_id = None
        job.status = "P"
        job.submit()

    def __resolve_job(self, job, queue_status):
        if job.status == "Q":
            logger.debug("%s is still waiting in queue.", job.name)
            self.current_jobs_statuses["Q"] += 1
        elif job.status == "R":
            logger.debug("%s is running with reported time: %s", job.name, queue_status[job.name]["time"])
            self.current_jobs_statuses["R"] += 1
        elif job.status == "D":
            logger.info("%s is no longer running verifying results", job.name)
            if job.is_completed():
                job.status = "S"
            else:
                self.num_failed_jobs += 1
                self.__resubmit_failed_job(job)
                self.current_jobs_statuses["Q"] += 1

        if job.status == "S":
            job.synchronize(self.datadir_local)
            job.status = "F"
            self.current_jobs_statuses["F"] += 1

    def __derive_jobs_info_from_files(self):
        for job in self.job_list:
            if job.status == "F":
                self.current_jobs_statuses["F"] += 1
            else:
                logger.info("%s status unknown verifying results", job.name)
                # it is either done or has crushed
                if job.is_completed():
                    job.synchronize(self.datadir_local)
                    job.status = "F"
                    self.current_jobs_statuses["F"] += 1
                else:
                    self.num_failed_jobs += 1
                    self.__resubmit_failed_job(job)
                    self.current_jobs_statuses["Q"] += 1

    def __update_template_files(self):
        file_template = create_docking_batch_script("vina_template.txt")
        logger.info("remote.update_template_file = True >>> updating batch scripts \"run_docking.sh\"")

        for folder in os.listdir(self.datadir_local):
            if "batch" not in folder:
                continue
            file_path = os.path.join(self.datadir_local, folder, "run_docking.sh")
            with open(file_path, "w") as output_stream:
                output_stream.write(file_template.replace("$$REMOTE.BATCH_NAME$$", folder))

    def __prepare_local_data(self):
        if not os.path.exists(self.datadir_local):
            raise RuntimeError("Local folder with batches {} does not exist."
                               "run script for candidates preparation to get it.".format(self.datadir_local))

        tmp_cwd = os.getcwd()
        os.chdir(self.datadir_local)
        logger.info("Creating local compressed file with data for transfer.")
        with tarfile.open(self.dataname_local + ".tar.gz", "w:gz") as transport_tar:
            for folder in os.listdir(self.datadir_local):
                if "batch" not in folder:
                    continue
                transport_tar.add(folder)
        os.chdir(tmp_cwd)

    def __transfer_data_to_remote(self):
        if self.__remote_datadir_exists():
            raise RuntimeError("Screening directory {} already exists on remote host, rename it or "
                               "alter the parameter \"docking_dataname\" in [REMOTE] section of "
                               "the configuration file".format(self.dataname_remote))
        else:
            self.remote_connection.makedir(self.remote_screening_dir)
            if not self.__remote_datadir_exists():
                raise RuntimeError("Failed to create screening directory {} "
                                   "on the remote host".format(self.dataname_remote))

        self.remote_connection.send_file(self.source, self.dest)
        file_local = os.stat(self.source)
        file_remote = self.remote_connection.file_info(self.dest)

        if (file_remote is not None) and (round(file_local.st_mtime, 0) <= round(file_remote.st_mtime, 0)):
            logger.info("Package with screening data successfully transferred to remote host. Unpacking")
            utils.remove_existing(self.source)
            commands = ["cd {} && tar xzf {}.tar.gz".format(self.remote_screening_dir, self.dataname_local)]
            self.remote_connection.execute_remote_command(commands)
            self.remote_connection.remove_file(self.dest)
        else:
            raise RuntimeError("Could not set up screening on the remote host.")

    def __set_jobs(self):
        for batch_dir in sorted(self.remote_connection.listdir(self.remote_screening_dir)):
            if "batch" in batch_dir:
                self.job_list.append(RemoteJob(os.path.join(self.remote_screening_dir, batch_dir),
                                               self.remote_connection, self.queue_manager, self.names_to_delete))

    def set_up_screening(self):
        if os.path.exists(self._screening_status_storage):
            logger.warning("Saved screening found, restarting from file: %s", self._screening_status_storage)
            logger.warning("Should you like to start new screening, please remove this file")
            self.__load_screening()
        else:
            logger.info("Setting the remote screening folder at %s@%s", REMOTE["host"]["username"],
                        REMOTE["host"]["hostname"])

            if utils.is_candidate_processing_status("docking-failed") and REMOTE["update_template_file"]:
                self.__update_template_files()
                utils.update_candidate_processing_status("docking-prepared")

            if utils.is_candidate_processing_status("docking-prepared"):
                self.__prepare_local_data()
                utils.update_candidate_processing_status("ready_for_transfer2remote")

            logger.info("Local data ready, proceeding with transfer")
            if utils.is_candidate_processing_status("ready_for_transfer2remote"):
                self.__transfer_data_to_remote()
                utils.update_candidate_processing_status("transferred2remote")

            if utils.is_candidate_processing_status("transferred2remote"):
                logger.info("Setting screening jobs from data on the remote host.")
                self.__set_jobs()
                self.__save_screening()

    def start_screening(self):
        if not (self.screening_finished or self.__is_frequent_job_error()):
            logger.info("Starting screening with %d batches", len(self.job_list))
            counter = 0
            for job in self.job_list:

                if job.job_id is not None:
                    logger.info("Reloading job with ID %s for %s", job.job_id, job.name)
                else:
                    job.submit()
                    time.sleep(0.5)
                    counter += 1
                    if counter % 100 == 0:
                        self.__save_screening()
                        logger.info("%d/%d jobs submitted", counter, len(self.job_list))

            self.__save_screening()

    def check_jobs(self):
        if self.screening_finished or self.__is_frequent_job_error():
            return True

        self.current_jobs_statuses = {
            "F": 0,
            "R": 0,
            "Q": 0
        }

        queue_status = self.__get_queue_status()

        if queue_status is None:
            logger.warning("Queue manager status inaccessible, will try again later")
            time.sleep(60)
            return False
        elif queue_status == ['']:
            logger.warning("Queue manager status does not contain info on batches")
            self.__derive_jobs_info_from_files()
            return self.__evaluate_screening_status()

        for job in self.job_list:
            job.sync_ids(queue_status)
            if job.status == "F":
                self.current_jobs_statuses["F"] += 1
            else:
                if job.job_id is not None:
                    if job.name in queue_status.keys():
                        job.status = queue_status[job.name]["status"]
                    else:
                        # job was submitted and we have no info on its status from the queue manager
                        job.status = "D"

                    self.__resolve_job(job, queue_status)
                else:
                    # no job_id >> job was never successfully submitted to queue
                    job.clean_job_dir()
                    job.status = "P"
                    job.submit()

        return self.__evaluate_screening_status()

    def sync_data_locally(self):
        self.failed_files = {}
        if not self.screening_finished:
            return

        verification_tasks = _create_docking_verification_tasks(self.datadir_local)

        pool = mp.Pool(processes=ENV["num_parallel_cpu"])
        processing = []

        for task in verification_tasks:
            processing.append(pool.apply_async(_verify_docking_result, args=task))

        for p in processing:
            batch_name, failed_files_in_batch = p.get()
            if len(failed_files_in_batch) != 0:
                self.failed_files[batch_name] = failed_files_in_batch
        pool.close()

        if self.failed_files:
            self.__analyze_failed_files()
        else:
            logger.info("Removing remote screening data")
            self.__delete_remote_screening_data()
            logger.info("Removing local temporary data")
            shutil.rmtree(self.datadir_local)
            utils.remove_existing(self._screening_status_storage)
            utils.update_candidate_processing_status("docking-completed")

    def __resolve_failed_files(self):
        resolved_files = {}
        failed_cids = []
        for batch_name in self.failed_files.keys():
            for failed_file in self.failed_files[batch_name]:
                failed_cids.append(failed_file[1][:-len(".pdbqt")])

        for batch_name in list(self.failed_files.keys()):
            path_to_dir = os.path.join(self.remote_screening_dir, batch_name)
            batch_output = {}
            for file_ in self.remote_connection.listdir(path_to_dir):
                if self.output_filename_pattern in file_:
                    file_path = os.path.join(self.remote_screening_dir, batch_name, file_)
                    previous_line = None
                    ligand_cid = None

                    for line in self.remote_connection.read_file(file_path):
                        if "#################################################################" in line:
                            if "#" not in previous_line and previous_line is not None:
                                past_ligand_cid = ligand_cid
                                ligand_cid = previous_line.rstrip()
                                if ligand_cid in failed_cids:
                                    if past_ligand_cid is not None:
                                        batch_output[past_ligand_cid]["output"].pop()

                                    batch_output[ligand_cid] = {"output": [previous_line],
                                                                "errors": False}

                        if ligand_cid is not None and ligand_cid in failed_cids:
                            batch_output[ligand_cid]["output"].append(line)
                            if "error" in line.lower() or "warn" in line.lower():
                                batch_output[ligand_cid]["errors"] = True

                        previous_line = line

            for failed_file in self.failed_files[batch_name].copy():
                subdir = failed_file[0]
                file_name = failed_file[1]
                cid = file_name[:-len(".pdbqt")]
                if cid in batch_output.keys() and batch_output[cid]["errors"]:
                    logger.info("Error resolved for failed file %s in %s.", file_name, subdir)
                    for line in batch_output[cid]["output"]:
                        logger.debug("The output follows:" + line)
                    if batch_name not in resolved_files.keys():
                        resolved_files[batch_name] = {}
                    if subdir not in resolved_files[batch_name].keys():
                        resolved_files[batch_name][subdir] = {}
                    if file_name not in resolved_files[batch_name][subdir].keys():
                        resolved_files[batch_name][subdir][file_name] = []
                        resolved_files[batch_name][subdir][file_name].append(failed_file)
                    if failed_file in self.failed_files[batch_name]:
                        self.failed_files[batch_name].remove(failed_file)

        if len(self.failed_files[batch_name]) == 0:
            del self.failed_files[batch_name]

        if resolved_files:
            output_file = os.path.join(paths.candidates(), "resolved_failed_candidates.json")
            with open(output_file, "w") as out_stream:
                json.dump(resolved_files, out_stream)

    def __clean_local_failed_batches(self):
        for batch in self.failed_files.keys():
            batch_path = os.path.join(self.datadir_local, batch)
            for folder in [os.path.join(batch_path, "log"), os.path.join(os.path.join(batch_path, "out"))]:
                for subdir in os.listdir(folder):
                    try:
                        os.rmdir(os.path.join(folder, subdir))
                    except OSError:
                        pass

    def __count_failed_files(self):
        num_failed_files = 0
        for batch_name in self.failed_files.keys():
            num_failed_files += len(self.failed_files[batch_name])
        return num_failed_files

    def __report_unresolved_files(self):
        output_file = os.path.join(paths.candidates(), "unresolved_failed_candidates.txt")
        with open(output_file, "w") as out_stream:
            for batch_name in self.failed_files.keys():
                for failed_file in self.failed_files[batch_name]:
                    out_stream.write("{}:: {}: {}\n".format(batch_name, failed_file[0], failed_file[1]))

    def __analyze_failed_files(self):
        if self.__remote_datadir_exists():
            num_failed_files = self.__count_failed_files()
            logger.warning("Docking of %d files failed, you might see the reasons bellow", num_failed_files)
            self.__resolve_failed_files()

            num_unresolved_failed_files = self.__count_failed_files()
            if num_unresolved_failed_files == 0:
                logger.info("Find at least one error msg for all failed files")
                logger.info("Removing remote screening data")
                self.__delete_remote_screening_data()
                logger.info("Removing local temporary data")
                shutil.rmtree(self.datadir_local)
                utils.remove_existing(self._screening_status_storage)
                utils.update_candidate_processing_status("docking-failed")
            else:
                logger.warning("Origin of failure remains unresolved for %d/%d files",
                               num_unresolved_failed_files, num_failed_files)
                self.__report_unresolved_files()
                logger.info("Removing resolved remote screening data")
                self.__delete_remote_screening_data()
                self.__clean_local_failed_batches()
                utils.update_candidate_processing_status("docking-failed")


class RemoteJob:
    def __init__(self, folder_path, handler, queue_manager, names_to_delete):
        self.location = folder_path
        self.queue_manager = queue_manager
        self.name = os.path.basename(folder_path)
        self.status = "P"
        # statuses: Prepared, Queued, Running, Done, Success, Finished
        self.job_id = None
        self.names_to_delete = names_to_delete
        self.remote_connection = handler

    def sync_ids(self, queue_status):
        if self.name in queue_status.keys() and self.job_id != queue_status[self.name]["job_id"]:
            logger.debug("For %s new job_ID found, updating %s >> %s",
                         self.name, self.job_id, queue_status[self.name]["job_id"])
            self.job_id = queue_status[self.name]["job_id"]

    def submit(self):
        if self.job_id is not None:
            logger.warning("Job in %s has already been submitted running with ID %s and has not been "
                           "further processed.", self.location, self.job_id)
        else:
            commands = ["cd {} && {} run_docking.sh".format(self.location, self.queue_manager.get_submit_job_command())]
            output = self.remote_connection.execute_remote_command(commands).decode('ascii').split("\n")
            self.job_id = output[0].split("job")[-1].strip()
            self.status = "Q"
            logger.info("Starting job with ID %s for %s", self.job_id, self.name)

    def report(self):
        return {
            "name": self.name,
            "status": self.status,
            "job_id": self.job_id,
            "location": self.location
        }

    def clean_job_dir(self):
        if self.remote_connection.file_info(self.location) is not None:
            paths2remove = []
            filenames_to_delete = self.names_to_delete.copy()
            for considered_file in self.remote_connection.listdir(self.location):
                for filename in filenames_to_delete:
                    if filename in considered_file:
                        filenames_to_delete.remove(filename)
                        paths2remove.append(os.path.join(self.location, considered_file))
                        break
            self.remote_connection.remove_files(paths2remove)

    def synchronize(self, datadir_local):
        logger.info("%s docked successfully, copying results to local directories", self.name)
        for filename in ["out.tar.gz", "log.tar.gz"]:
            source = os.path.join(self.location, filename)
            dest = os.path.join(datadir_local, self.name, filename)
            self.remote_connection.get_file(source, dest)
            if filename.endswith(".tar.gz"):
                with tarfile.open(dest, "r:gz") as tar:
                    tar.extractall(os.path.join(datadir_local, self.name))
                utils.remove_existing(dest)

    def is_completed(self):
        remote_out_tarfile = self.remote_connection.file_info(os.path.join(self.location, "out.tar.gz"))
        remote_log_tarfile = self.remote_connection.file_info(os.path.join(self.location, "log.tar.gz"))
        if (remote_out_tarfile is not None) and (remote_log_tarfile is not None):
            if (remote_log_tarfile.st_size > 0) and (remote_out_tarfile.st_size > 0):
                return True

        return False


def _create_docking_verification_tasks(datadir_local):
    tasks = []
    for batch_name in os.listdir(datadir_local):
        batch_dir = os.path.join(datadir_local, batch_name)
        with tarfile.open(os.path.join(batch_dir, "pdbqt.tar.gz"), "r:gz") as tar:
            file_names = tar.getnames()
        tasks.append([{"file_names": file_names,
                       "batch_dir": batch_dir,
                       "batch_name": batch_name,
                       "datadir_local": datadir_local
                       }])
    return tasks


def _verify_docking_result(task):
    batch_name = task["batch_name"]
    batch_dir = task["batch_dir"]
    failed_files_in_batch = []
    logger.info("Verifying results in batch: %s", batch_name)
    for name in task["file_names"]:
        subdir, pdbqt_file = paths.get_subdir_filename_from_filepath(name)
        logger.debug("Testing results of file %s in %s/%s", pdbqt_file, batch_name, subdir)

        outfile = os.path.join(batch_dir, "out", subdir, pdbqt_file)
        logfile = os.path.join(batch_dir, "log", subdir, pdbqt_file.replace(".pdbqt", ".log"))
        dest_outfile = os.path.join(paths.candidates_docking(), subdir, pdbqt_file)
        dest_logfile = os.path.join(paths.candidates_docking(), subdir, pdbqt_file.replace(".pdbqt", ".log"))

        if os.path.exists(dest_outfile) and os.path.exists(dest_logfile):
            utils.remove_existing(outfile)
            utils.remove_existing(logfile)
        else:
            if test_vina_results(logfile, outfile):
                logger.debug("Consistent results found for file %s in %s ", pdbqt_file, subdir)
                os.makedirs(os.path.join(paths.candidates_docking(), subdir), exist_ok=True)
                shutil.move(outfile, dest_outfile)
                shutil.move(logfile, dest_logfile)
            else:
                logger.warning("Could not verify results for file %s in %s", pdbqt_file, subdir)
                failed_files_in_batch.append((subdir, pdbqt_file))

    if len(failed_files_in_batch) == 0:
        logger.info("%s was successful completely, removing its temporary directory from %s.",
                    batch_name, task["datadir_local"])
        shutil.rmtree(batch_dir)

    return batch_name, failed_files_in_batch


