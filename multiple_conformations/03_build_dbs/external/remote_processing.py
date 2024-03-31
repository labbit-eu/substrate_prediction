#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import paramiko
import logging
import socket
import time
import stat

from libs.configuration import REMOTE

logger = logging.getLogger(__name__)


class RemoteHandler:
    def __init__(self):
        self.ssh_output = None
        self.ssh_error = None
        self.client = None
        self.host = REMOTE["host"]["hostname"]
        self.username = REMOTE["host"]["username"]
        self.port = int(REMOTE["host"]["port"])
        self.ssh_key = paramiko.RSAKey.from_private_key_file(REMOTE["host"]["ssh_key_path"])
        logging.getLogger("paramiko").setLevel(logging.WARN)

    def __connect(self):
        try:
            self.client = paramiko.SSHClient()
            self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            self.client.connect(hostname=self.host, username=self.username, pkey=self.ssh_key)
        except paramiko.AuthenticationException:
            raise RuntimeError("Authentication failed, please verify your credentials")
        except paramiko.SSHException:
            raise RuntimeError("Could not establish SSH connection: ", paramiko.SSHException)
        except socket.timeout as err:
            raise RuntimeError("Connection timed out: ", err)
        except Exception as exc:
            raise RuntimeError("Exception in connecting to the server: ", exc)

    def makedir(self, remote_path):
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()
            results = sftp_client.mkdir(remote_path)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()
        return results

    def listdir(self, remote_path):
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()
            results = sftp_client.listdir(remote_path)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()
        return results

    def make_folder_path(self, path):
        status_path = self.file_info(path)
        if status_path is not None:
            if stat.S_ISDIR(status_path.st_mode):
                return True
            else:
                return False
        else:
            path, folder = os.path.split(path)
            root_dir_exists = self.make_folder_path(path)
            if root_dir_exists:
                self.makedir(os.path.join(path, folder))
                return True
            else:
                return False

    def remove_jobdir(self, path):
        status_path = self.file_info(path)
        if status_path is not None and stat.S_ISDIR(status_path.st_mode):
            try:
                self.__connect()
                sftp_client = self.client.open_sftp()
            except Exception as exc:
                raise RuntimeError("Following exception reported: ", exc)

            for folder in ["out", "log"]:
                folder_path = os.path.join(path, folder)
                folder_status = self.file_info(folder_path)
                if folder_status is not None and stat.S_ISDIR(folder_status.st_mode):
                    for subdir in sftp_client.listdir(folder_path):
                        subdir_path = os.path.join(folder_path, subdir)
                        sftp_client.rmdir(subdir_path)
                    sftp_client.rmdir(folder_path)

            folder_path = os.path.join(path, "in")
            folder_status = self.file_info(folder_path)
            if folder_status is not None and stat.S_ISDIR(folder_status.st_mode):
                for file_name in sftp_client.listdir(folder_path):
                    file_path = os.path.join(folder_path, file_name)
                    sftp_client.remove(file_path)

                sftp_client.rmdir(folder_path)

            sftp_client.rmdir(path)
            sftp_client.close()
            self.client.close()

    def remove_files(self, paths):
        if paths:
            try:
                self.__connect()
                sftp_client = self.client.open_sftp()
            except Exception as exc:
                raise RuntimeError("Following exception reported: ", exc)

            for path in paths:
                sftp_client.remove(path)

            sftp_client.close()
            self.client.close()

    def remove_dir(self, remote_path):
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()

            sftp_client.rmdir(remote_path)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()

    def remove_file(self, remote_path):
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()
            sftp_client.remove(remote_path)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()

    def file_info(self, remote_path):
        results = None
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()
            if os.path.basename(remote_path) in sftp_client.listdir(os.path.dirname(remote_path)):
                results = sftp_client.stat(remote_path)
            else:
                logger.warning("file %s does not exist at remote location.", remote_path)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()
            return results

    def read_file(self, remote_path):
        results = None
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()
            if os.path.basename(remote_path) in sftp_client.listdir(os.path.dirname(remote_path)):
                remote_stream = sftp_client.open(remote_path)
                results = remote_stream.read().decode('ascii').split("\n")
                remote_stream.close()
            else:
                logger.warning("file %s does not exist at remote location.", remote_path)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()
            return results

    def send_file(self, local_path, remote_path, mode=None):
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()
            sftp_client.put(local_path, remote_path)
            if mode is not None:
                sftp_client.chmod(remote_path, mode)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()

    def get_file(self, remote_path, local_path):
        try:
            self.__connect()
            sftp_client = self.client.open_sftp()
            sftp_client.get(remote_path, local_path)
        except Exception as exc:
            raise RuntimeError("Following exception reported: ", exc)
        finally:
            sftp_client.close()
            self.client.close()

    def execute_remote_command(self, commands, repeat=False):
        self.ssh_output = None
        try:
            self.__connect()
            for command in commands:
                logger.debug("Executing command --> {}".format(command))
                stdin, stdout, stderr = self.client.exec_command(command, timeout=3600)
                self.ssh_output = stdout.read()
                self.ssh_error = stderr.read()
                if self.ssh_error:
                    if repeat:
                        count_errors = 1
                        time.sleep(30)
                        while self.ssh_error and count_errors <= 5:
                            logger.error("Error %s occurred while running command: %s", self.ssh_error.decode('ascii'),
                                         command)
                            time.sleep(count_errors * 30)
                            stdin, stdout, stderr = self.client.exec_command(command, timeout=3600)
                            self.ssh_output = stdout.read()
                            self.ssh_error = stderr.read()
                            count_errors += 1

                        if self.ssh_error:
                            raise RuntimeError("Persistent error occurred %d-times".format(count_errors))
                    else:
                        logger.error("Error %s occurred while running command: %s", self.ssh_error.decode('ascii'),
                                     command)
        except socket.timeout:
            raise RuntimeError("Command(s) %s timed out.".format(commands))
        except paramiko.SSHException:
            raise RuntimeError("Failed to execute command(s): ", commands)
        finally:
            self.client.close()
            return self.ssh_output


#  queue managers
# TODO test if the queueManager is of supported version!!
class QueueManager:
    def __init__(self):
        self.type = REMOTE["queue"]["manager"]
        self.status_converter = None
        self.status_command = None
        self.submit_command = None
        self.output_filename_pattern = None

    def get_queue_status_command(self):
        return self.status_command

    def get_submit_job_command(self):
        return self.submit_command

    def get_output_filename_pattern(self):
        return self.output_filename_pattern

    def get_handler(self):
        if self.type == "slurm":
            return _SlurmHandler()
        elif self.type == "torque":
            return _TorqueHandler()


class _SlurmHandler(QueueManager):
    def __init__(self):
        self.status_converter = {
            "BF": "D",
            "CA": "D",
            "CD": "D",
            "CF": "Q",
            "CG": "R",
            "DL": "D",
            "F": "D",
            "NF": "D",
            "OOM": "D",
            "PD": "Q",
            "PR": "D",
            "R": "R",
            "RD": "Q",
            "RH": "Q",
            "RQ": "Q",
            "SE": "Q",
            "ST": "D",
            "S": "Q",
            "TO": "D"
        }
        self.status_command = "squeue -h --format=\"%.18i %.9P %.10j %.8u %.2t %.10M %.6D %R\" -u {}".format(
            REMOTE["host"]["username"])
        self.submit_command = "sbatch"
        self.output_filename_pattern = "slurm"

    def parse_job_info_line(self, line):
        remote_status = line.split()[4]
        if remote_status not in self.status_converter.keys():
            raise KeyError("Reported job status %s unknown", remote_status)

        return {
            "job_id": line.split()[0],
            "queue": line.split()[1],
            "name": line.split()[2],
            "user": line.split()[3],
            "status": self.status_converter[remote_status],
            "time": line.split()[5]
            }


class _TorqueHandler(QueueManager):
    def __init__(self):
        self.status_converter = {
            "C": "D",
            "E": "R",
            "H": "Q",
            "Q": "Q",
            "R": "R",
            "T": "Q",
            "W": "Q",
        }
        self.status_command = "qstat -u {}".format(REMOTE["host"]["username"])
        self.submit_command = "qsub"
        self.output_filename_pattern = ""

    def parse_job_info_line(self, line):
        remote_status = line.split()[4]
        if remote_status not in self.status_converter.keys():
            raise KeyError("Reported job status %s unknown", remote_status)

        return {
            "job_id": line.split()[0],
            "queue": line.split()[5],
            "name": line.split()[1],
            "user": line.split()[2],
            "status": self.status_converter[remote_status],
            "time": line.split()[3]
            }
