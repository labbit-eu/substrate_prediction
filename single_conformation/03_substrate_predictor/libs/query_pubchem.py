#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Properties:
 * HBondAcceptor_value
 * HBondDonor_value
 * HeavyAtomCount_value
 * MolecularWeight_value
 * RotatableBond_value
 * XLogP3_value
 * canSMILES_value
 * isoSMILES_value
"""

import math
import os
import time
import re
import io
import logging
import gzip
import abc
import tempfile
import json
import ftplib
import itertools
import collections
import multiprocessing
import requests
from datetime import datetime
import urllib.request

from libs.utils import remove_existing, read_json
from libs import paths
from libs.configuration import ENV, CHEM, TASK

from rdkit import Chem

FILTER_MERGE_STEP_SIZE = 5000000
SUBSET_SIZE = 10000

logger = logging.getLogger(__name__)


class Ftp(object):
    """Mirror remote files from FTP to local host."""

    def __init__(self, host, username="anonymous", password="anonymous"):
        self._host = host
        self._username = username
        self._password = password

    def mirror_directory(self, remote, pattern, local):
        os.makedirs(local, exist_ok=True)
        pattern = re.compile(pattern)
        with ftplib.FTP(self._host, self._username, self._password) as ftp:
            ftp.cwd(remote)
            logger.debug("Listing files in: %s", remote)
            files = [file for file in ftp.nlst()
                     if pattern.match(file)]
            for index, file in enumerate(files):
                output_path = os.path.join(local, file)
                if os.path.exists(output_path):
                    logger.debug("Skipping %d/%d %s -> %s",
                                 index + 1, len(files), file, output_path)
                    continue
                logger.debug("Copying %d/%d %s -> %s",
                             index + 1, len(files), file, output_path)
                self._retrieve_file(ftp, file, output_path)

    def mirror_file(self, remote, local):
        logger.debug("Copying %s -> %s", remote, local)
        with ftplib.FTP(self._host, self._username, self._password) as ftp:
            remote_dir = remote[:remote.rindex("/")]
            remote_file = remote[remote.rindex("/") + 1:]
            ftp.cwd(remote_dir)
            self._retrieve_file(ftp, remote_file, local)

    @staticmethod
    def _retrieve_file(ftp, remote, local):
        command = "RETR " + remote
        with open(local, "wb") as stream:
            ftp.retrbinary(command, stream.write)


class Sort(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def consume(self, key, value):
        raise NotImplementedError()

    @abc.abstractmethod
    def sort(self):
        raise NotImplementedError()

    @abc.abstractmethod
    def iterate(self):
        raise NotImplementedError()


class InMemorySort(Sort):

    def __init__(self):
        self._values = []

    def consume(self, key, value):
        self._values.append([key, value])

    def sort(self):
        self._values.sort(key=lambda tuple: tuple[0])

    def iterate(self):
        for key, value in self._values:
            yield key, value


class ExternalMergeSort(Sort):

    def __init__(self, bucket_size=5000000, parallel_files=64):
        self._parallel_files = parallel_files
        self._bucket_size = bucket_size
        #
        self._values = []
        self._directory = tempfile.mkdtemp(None, "python_merge_sort-")
        os.makedirs(self._directory, exist_ok=True)
        self._files = []
        self._files_counter = 0

    def consume(self, key, value):
        self._values.append((key, value))
        if len(self._values) >= self._bucket_size:
            self._flush()

    def _flush(self):
        self._values.sort(key=lambda tuple: tuple[0])
        file = self._new_file()
        with gzip.open(file, "wb") as stream:
            for key, value in self._values:
                stream.write("{},{}\n".format(key, value).encode("utf-8"))
        self._values = []
        self._files.append(file)

    def _new_file(self):
        self._files_counter += 1
        return os.path.join(
            self._directory,
            "{}.txt.gz".format(self._files_counter))

    def sort(self):
        self._flush()
        while len(self._files) > 1:
            logger.debug("Reducing files: %d", len(self._files))
            split_index = min(self._parallel_files, len(self._files))
            files_to_merge = self._files[0:split_index]
            self._files = self._files[split_index:]
            self._files.append(self._reduce(files_to_merge))

    def _reduce(self, files):
        path = self._new_file()
        sources = []
        for file in files:
            sources.append(self._consume_file(file))
        values = [next(iterator) for iterator in sources]
        with gzip.open(path, "wb") as stream:
            while len(values) > 0:
                # TODO Optimize minimum selection.
                value = min(values, key=lambda tuple: int(tuple[0]))
                index = values.index(value)
                try:
                    values[index] = next(sources[index])
                except StopIteration:
                    del values[index]
                    del sources[index]
                stream.write("{},{}\n".format(*value).encode("utf-8"))
        return path

    def _consume_file(self, file):
        with gzip.open(file, "rb") as stream:
            for line in stream:
                line = line.decode("utf-8")
                key, value = line.rstrip().split(",", maxsplit=2)
                yield int(key), value
        remove_existing(file)

    def iterate(self):
        assert len(self._files) == 1, "Data were not sorted."
        yield from self._consume_file(self._files[0])
        os.removedirs(self._directory)


class Executor(object):
    """Execute function with given arguments. Supports multiprocessing."""

    def __init__(self, multiprocess, processes):
        self._processes = processes
        self._pool = None
        self._multiprocess = multiprocess
        self._function = None

    def execute(self, funct, arguments):
        if self._multiprocess:
            self.execute_multiprocess(funct, arguments)
        else:
            self.execute_single_process(funct, arguments)

    @staticmethod
    def execute_single_process(funct, arguments):
        logger.info("Executing %d tasks in a single tread.", len(arguments))
        for args in arguments:
            funct(args)

    def execute_multiprocess(self, funct, arguments):
        logger.info("Executing %d tasks using %d processes.",
                    len(arguments), self._processes)
        if self._pool is None:
            self._pool = multiprocessing.Pool(processes=self._processes)
        self._pool.map(funct, arguments)

    def close(self):
        if self._pool is not None:
            self._pool.close()
            self._pool = None

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    @staticmethod
    def sequential():
        return Executor(multiprocess=False, processes=None)

    @staticmethod
    def parallel(processes_count):
        return Executor(multiprocess=True, processes=processes_count)


class PubchemTurtleFile(object):

    def __init__(self, path):
        self._path = path
        self._stream = None
        self._reader = None

    def __enter__(self):
        self._stream = gzip.open(self._path, "rb")
        self._reader = io.BufferedReader(self._stream)
        return self

    def __exit__(self, type, value, traceback):
        if self._reader is not None:
            self._reader.close()
            self._reader = None
        if self._stream is not None:
            self._stream.close()
            self._stream = None

    def __iter__(self):
        return self

    def __next__(self):
        for line in self._reader:
            if b"has-value" not in line:
                continue
            line = line.decode("utf-8")
            cid = line[line.find(":") + 4: line.find("_")]
            value = line[line.rfind("\t") + 2: line.rfind("\"")]
            return int(cid), value
        raise StopIteration


class WritableValueFile(object):

    def __init__(self, path, info=None):
        self._path = path
        self._info = info
        self._stream = None
        self._writer = None
        self._min = None
        self._max = None
        self._count = 0

    def __enter__(self):
        self._stream = gzip.open(self._path, "wt")
        return self

    def __exit__(self, type, value, traceback):
        if self._stream is not None:
            self._stream.close()
            self._stream = None
        with open(self._path + ".json", "w") as stream:
            stats = {
                "min": self._min,
                "max": self._max,
                "count": self._count
            }
            if self._info is not None:
                stats.update(self._info)
            json.dump(stats, stream)

    def consume(self, key, value):
        self._count += 1
        if self._min is None or self._min > key:
            self._min = key
        if self._max is None or self._max < key:
            self._max = key
        self._stream.write("{},{}\n".format(key, value))


class ReadableValueFile(object):

    def __init__(self, path):
        self._path = path
        self._stream = None
        self._reader = None

    def __enter__(self):
        self._stream = gzip.open(self._path, "rt")
        return self

    def __exit__(self, type, value, traceback):
        if self._stream is not None:
            self._stream.close()
            self._stream = None

    def __iter__(self):
        return self

    def __next__(self):
        for line in self._stream:
            line = line.rstrip().split(",", maxsplit=2)
            return int(line[0]), line[1]
        raise StopIteration


def query_pubchem_for_substrates(filters):
    rdf_dir = os.path.join(".", "cache", "pubchem", "rdf")
    _download_files_with_values(filters.keys(), rdf_dir)

    prop_dir = paths.candidates_pubchem()
    temp_dir = os.path.join(paths.temp(), "pubchem_properties")

    with Executor.parallel(ENV["num_parallel_cpu"]) as executor:
        logger.info("Creating filtering tasks ...")
        filter_tasks = _create_filter_rdf_tasks(filters, rdf_dir, prop_dir)

        logger.info("Filtering ...")
        executor.execute(_filter_rdf_files, filter_tasks)
        filter_tasks += _get_cids_with_structure(prop_dir, executor)

        logger.info("Creating merging tasks ...")
        merger_tasks = _create_tasks(prop_dir, filter_tasks, FILTER_MERGE_STEP_SIZE)

        logger.info("Merging ...")
        executor.execute(_merge_value_files, merger_tasks)

        number_candidates_for_download, cids_to_get = _count_candidates(prop_dir)
        if number_candidates_for_download > TASK["candidates"]["prepare"]["max_molecules"]:
            logger.warning("Too many candidates (%d) for further processing. Reducing their number to fit %d set "
                           "in the parameter \"candidates.prepare.max_molecules\" in [TASK] section of "
                           "the configuration file.", number_candidates_for_download,
                           TASK["candidates"]["prepare"]["max_molecules"])
            cids_to_get = _prune_candidates(filters, executor, prop_dir, temp_dir)

        number_candidates_for_download = _update_input_files_for_download(prop_dir, paths.candidates_sdf(), cids_to_get)

        logger.info("Downloading %d candidates ...", number_candidates_for_download)
        download_tasks = _create_download_tasks(prop_dir, paths.candidates())

        executor.execute(_download_from_pubchem_by_cid, download_tasks)
        downloaded_candidates = _evaluate_sdf_output_dir(paths.candidates_sdf())

        if downloaded_candidates < number_candidates_for_download:
            logger.warning("%d candidates have not been downloaded. Trying once more.",
                           number_candidates_for_download - downloaded_candidates)
            executor.execute(_download_from_pubchem_by_cid, download_tasks)
            downloaded_candidates = _evaluate_sdf_output_dir(paths.candidates_sdf())

    if downloaded_candidates < number_candidates_for_download:
        logger.error("%d/%d candidates have not been downloaded from PubChem.",
                     number_candidates_for_download - downloaded_candidates, number_candidates_for_download)
    else:
        logger.info("Done... %d/%d candidates have been downloaded.", downloaded_candidates,
                    number_candidates_for_download)


def _prune_candidates(filters, executor, prop_dir, temp_dir):
    if os.path.exists(os.path.join(temp_dir, "property_file_final.txt.gz")) \
            and os.path.exists(os.path.join(temp_dir, "distributions_final.json")):
        logger.debug("File with pruned candidates %s exists, loading.",
                     os.path.join(temp_dir, "property_file_final.txt.gz"))
        distributions = read_json(os.path.join(temp_dir, "distributions_final.json"))
        return _cids_to_get(distributions)

    tasks = _prepare_tasks_for_common_property_filtering(filters, prop_dir, temp_dir)
    executor.execute(_filter_property_files_for_common, tasks)
    logger.info("Merging property files")
    keys = _merge_common_to_single_file(temp_dir)
    logger.info("Merging done.. ")
    num_candidates = tasks[-1]["count"]

    distributions = _create_distributions(keys, temp_dir)
    final_distributions = _greedy_balanced_removal(distributions, num_candidates, temp_dir)
    return _cids_to_get(final_distributions)


def _get_e_utils_response_values(text, start_tag, end_tag):
    start_index = text.find(start_tag)
    end_index = text.find(end_tag)
    if end_index >= 0 and start_index >= 0:
        return text[start_index + len(start_tag):end_index]
    else:
        return False


def _prepare_tasks_for_get_cids_with_structure(prop_dir):
    count = None
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    query = "has_3d_conformer[filter]+AND+{}:{}[MW]+AND+{}:{}[RBC]+AND+{}:{}[XLGP]".format(
        TASK["candidates"]["prepare"]["MW"][0], TASK["candidates"]["prepare"]["MW"][1],
        TASK["candidates"]["prepare"]["TORSDOF"][0], TASK["candidates"]["prepare"]["TORSDOF"][1],
        TASK["candidates"]["prepare"]["logP"][0], TASK["candidates"]["prepare"]["logP"][1])
    url = base_url + "esearch.fcgi?db=pccompound&term={}&usehistory=y".format(query)
    while count is None:
        with urllib.request.urlopen(url) as response:
            body = response.read().decode('utf-8')
        key = _get_e_utils_response_values(body, "<QueryKey>", "</QueryKey>")
        web = _get_e_utils_response_values(body, "<WebEnv>", "</WebEnv>")
        count = int(_get_e_utils_response_values(body, "<Count>", "</Count>"))

    logger.info("Getting IDs of %d candidates with structures and %s between [%d,%d], %s between [%d,%d] and"
                " %s between [%d,%d]", count,
                "MolecularWeight", TASK["candidates"]["prepare"]["MW"][0], TASK["candidates"]["prepare"]["MW"][1],
                "RotatableBond", TASK["candidates"]["prepare"]["TORSDOF"][0],
                TASK["candidates"]["prepare"]["TORSDOF"][1],
                "XLogP3", TASK["candidates"]["prepare"]["logP"][0], TASK["candidates"]["prepare"]["logP"][1])

    retmax = 100000
    tasks = []
    for retstart in range(0, (math.ceil(count / retmax) + 1) * retmax, retmax):
        if retstart <= count:
            output_path = os.path.join(prop_dir, "cids_with_3D_structure_retstart_{}.txt.gz".format(retstart))
            task = {
                "base_url": base_url,
                "WebEnv": web,
                "QueryKey": key,
                "retmax": retmax,
                "retstart": retstart,
                "output_path": output_path
            }
            tasks.append(task)
    return tasks


def _download_cids_with_structure(task):
    output_path = task["output_path"]
    if os.path.exists(output_path + ".json"):
        logger.debug("File %s exists, skipping download", output_path + ".json")
        return

    min_cid = None
    max_cid = None
    counter = 0
    url = task["base_url"] + "esearch.fcgi?db=pccompound&WebEnv={}&query_key={}&retstart={}&retmax={}".format(
        task["WebEnv"], task["QueryKey"], task["retstart"], task["retmax"])

    time_windows = _assign_query_time(entrez=True)
    while not datetime.now().second in time_windows:
        time.sleep(0.6)

    logger.debug("Downloading cids with structures to %s", output_path)
    with urllib.request.urlopen(url) as response:
        body = response.read().decode('utf-8')
    ids = _get_e_utils_response_values(body, "<IdList>", "</IdList>")
    if ids:
        with gzip.open(output_path, "wt") as out_stream:
            for id_ in ids.split("\n"):
                cid = _get_e_utils_response_values(id_, "<Id>", "</Id>")
                if cid:
                    cid = int(cid)
                    if min_cid is None or cid < min_cid:
                        min_cid = cid
                    if max_cid is None or cid > min_cid:
                        max_cid = cid
                    counter += 1
                    out_stream.write("{},{}\n".format(cid, True))

        with open(output_path + ".json", "w") as json_stream:
            stats = {
                "min": min_cid,
                "max": max_cid,
                "count": counter
            }
            json.dump(stats, json_stream)


def _get_cids_with_structure(prop_dir, executor):
    cids_tasks = _prepare_tasks_for_get_cids_with_structure(prop_dir)
    executor.execute(_download_cids_with_structure, cids_tasks)
    tasks = []

    for file_name in os.listdir(prop_dir):
        if file_name.startswith("cids_with_3D_structure_retstart") and file_name.endswith(".txt.gz"):
            output_path = os.path.join(prop_dir, file_name)
            tasks.append({
                "property": "has_structure",
                "output": output_path
            })

    return tasks


def _prepare_tasks_for_common_property_filtering(filters, prop_dir, temp_dir):
    os.makedirs(temp_dir, exist_ok=True)
    candidates = {}

    for filter_file in sorted(os.listdir(prop_dir)):
        if "filter" not in filter_file or not filter_file.endswith(".txt.gz") or "subdir" in filter_file:
            continue
        filter_path = os.path.join(prop_dir, filter_file)
        with gzip.open(filter_path, "rt") as in_stream:
            for line in in_stream:
                cid = line.rstrip()
                subset_id = str(int(cid) // SUBSET_SIZE)
                if subset_id not in candidates.keys():
                    candidates[subset_id] = set()
                candidates[subset_id].add(cid)

    num_candidates = 0
    for i in candidates.keys():
        num_candidates += len(candidates[i])

    tasks = []
    for key in filters.keys():
        out_path = os.path.join(temp_dir, key + "_common.txt.gz")
        task = {
            "candidates": candidates,
            "count": num_candidates,
            "in_dir": prop_dir,
            "property": key,
            "output": out_path
        }
        tasks.append(task)
    return tasks


def _filter_property_files_for_common(task):
    logger.info("Enumerating %s of selected candidates", task["property"])
    if os.path.exists(task["output"]):
        logger.debug("Output file %s exists, skipping search.", task["output"])
    else:
        counter = 0
        with gzip.open(task["output"], "wt") as out_stream:
            for file_name in sorted(os.listdir(task["in_dir"])):
                if task["property"] not in file_name or not file_name.endswith(".txt.gz"):
                    continue
                file_path = os.path.join(task["in_dir"], file_name)
                logger.debug("reading file: %s", file_name)
                with gzip.open(file_path, "rt") as stream:
                    for line in stream:
                        cid, value = line.rstrip().split(",")
                        subset_id = str(int(cid) // SUBSET_SIZE)
                        if subset_id in task["candidates"].keys() and cid in task["candidates"][subset_id]:
                            counter += 1
                            out_stream.write("{},{}\n".format(cid, value))

        logger.info("Processed %d/%d molecules for %s", counter, task["count"], task["property"])


def _merge_common_to_single_file(temp_dir):
    out_file = os.path.join(temp_dir, "property_file_common.txt.gz")

    in_streams = {}
    keys = []

    for file_name in sorted(os.listdir(temp_dir)):
        if not file_name.endswith("_common.txt.gz") or "property_file_common.txt.gz" in file_name:
            continue
        key = file_name[:-len("_common.txt.gz")]
        keys.append(key)
        in_path = os.path.join(temp_dir, file_name)
        in_streams[key] = gzip.open(in_path, "rt")

    if os.path.exists(out_file):
        return keys

    with gzip.open(out_file, "wt") as out_stream:
        for line in in_streams[keys[0]]:
            cid, values = line.rstrip().split(",")
            for key in keys[1:]:
                line2 = in_streams[key].readline()
                values += "," + line2.rstrip().split(",")[1]
            out_stream.write("{},{}\n".format(cid, values))

    for key in keys:
        in_streams[key].close()

    return keys


def _create_distributions(keys, temp_dir):
    in_path = os.path.join(temp_dir, "property_file_common.txt.gz")
    storage = {}
    for key in keys:
        if "SMILES" in key:
            keys.remove(key)
        else:
            storage[key] = {"distribution": {}}

    out_path = os.path.join(temp_dir, "distributions_common.json")
    if os.path.exists(out_path):
        return read_json(out_path)

    with gzip.open(in_path, "rt") as in_stream:
        for line in in_stream:
            values = line.rstrip().split(",")[1:]
            for pos, key in enumerate(keys):
                value = str(int(float(values[pos])))
                if value not in storage[key]["distribution"].keys():
                    storage[key]["distribution"][value] = 0
                storage[key]["distribution"][value] += 1

    for key in keys:
        logger.info("Creating distribution of values for %s", key)
        values = [int(i) for i in storage[key]["distribution"].keys()]
        max_val = max(values)
        min_val = min(values)
        storage[key].update(
            {
                "min_val": min_val,
                "max_val": max_val,
                "pruning_step": math.ceil((max_val - min_val) / 20),
                "pruning_depth": 1
            }
        )
        storage.update(
            {
                "infile_path": in_path,
                "outfile_path": in_path.replace("_common.", "_pruned."),
                "properties": keys
            }
        )

    with open(out_path, "w") as out_stream:
        json.dump(storage, out_stream)

    return storage


def _update_distributions_source_file(distributions):
    in_path = distributions["infile_path"]
    out_path = distributions["outfile_path"]
    logger.debug("Pruning values in the distributions' source file %s", in_path)

    counter = 0
    storage = {}

    if "pruned" in in_path:
        out_path = distributions["outfile_path"] + ".bac"

    with gzip.open(out_path, "wt") as out_stream:
        with gzip.open(in_path, "rt") as in_stream:
            for line in in_stream:
                temp = line.rstrip().split(",")
                cid = temp[0]
                values = temp[1:-1]
                smiles = temp[-1]
                skip_line = False
                for key, value in zip(distributions["properties"], values):
                    storage[key] = {}
                    keep_values = distributions[key]["distribution"].keys()
                    value = str(int(float(value)))
                    if value not in keep_values:
                        skip_line = True
                if skip_line:
                    continue
                counter += 1
                out_values = ""
                for value in values:
                    out_values += value + ","

                out_stream.write("{},{}{}\n".format(cid, out_values, smiles))

    logger.debug("Updating the distributions... ")
    with gzip.open(out_path, "rt") as in_stream:
        for line in in_stream:
            values = line.rstrip().split(",")[1:]
            for pos, key in enumerate(distributions["properties"]):
                value = str(int(float(values[pos])))
                if value not in storage[key].keys():
                    storage[key][value] = 0
                storage[key][value] += 1

    for key in distributions["properties"]:
        distributions[key]["distribution"] = storage[key]

    if "pruned" in in_path:
        os.rename(out_path, in_path)
    else:
        distributions["infile_path"] = distributions["outfile_path"]
    logger.info("Keeping %d candidates... ", counter)
    return counter


def _cids_to_get(distributions):
    in_path = distributions["outfile_path"]
    cids_to_keep = {}
    with gzip.open(in_path, "rt") as in_stream:
        for line in in_stream:
            cid = line.rstrip().split(",")[0]
            subset_id = str(int(cid) // SUBSET_SIZE)
            if subset_id not in cids_to_keep.keys():
                cids_to_keep[subset_id] = []
            cids_to_keep[subset_id].append(cid)

    logger.info("Pruning done, keeping %d candidates with following properties:", distributions["count"])
    for key in distributions["properties"]:
        logger.info("%s: [%d, %d]", key, distributions[key]["min_val"], distributions[key]["max_val"])

    return cids_to_keep


def _remove_remaining_with_largest_values(distributions, to_remove, temp_dir):
    in_path = distributions["infile_path"]
    out_path = distributions["outfile_path"] + ".bac"

    cids_to_remove = set()
    properties = {}
    with gzip.open(in_path, "rt") as in_stream:
        for line in in_stream:
            temp = line.rstrip().split(",")
            cid = temp[0]
            values = temp[1:-1]
            for key, value in zip(distributions["properties"], values):
                if key not in properties.keys():
                    properties[key] = {}
                properties[key][cid] = float(value)

    tmp_list = []
    for key in properties.keys():
        for cid in sorted(properties[key], key=properties[key].__getitem__):
            tmp_list.append(cid)

    while len(cids_to_remove) < to_remove:
        cids_to_remove.add(tmp_list.pop())

    fast_cids_to_remove = {}
    for cid in cids_to_remove:
        subset_id = str(int(cid) // SUBSET_SIZE)
        if subset_id not in fast_cids_to_remove.keys():
            fast_cids_to_remove[subset_id] = []
        fast_cids_to_remove[subset_id].append(cid)

    counter = 0
    with gzip.open(out_path, "wt") as out_stream:
        with gzip.open(in_path, "rt") as in_stream:
            for line in in_stream:
                cid = line.rstrip().split(",")[0]
                subset_id = str(int(cid) // SUBSET_SIZE)
                if subset_id not in fast_cids_to_remove.keys() or cid not in fast_cids_to_remove[subset_id]:
                    counter += 1
                    out_stream.write(line)

    final_file = os.path.join(temp_dir, "property_file_final.txt.gz")
    os.rename(out_path, final_file)
    final_distributions = {"outfile_path": final_file,
                           "count": counter,
                           "properties": distributions["properties"]
                           }

    for key in properties.keys():
        final_distributions[key] = {}
        final_values = sorted(properties[key].values())
        final_distributions[key]["min_val"] = min(final_values)
        final_distributions[key]["max_val"] = max(final_values)

    dist_path = os.path.join(temp_dir, "distributions_final.json")
    with open(dist_path, "w") as out_stream:
        json.dump(final_distributions, out_stream)

    return final_distributions


def _select_property_values_for_removal(distributions, total_removed, num_to_remove):
    max_cost = {
        "key": None,
        "num_removed": 0,
        "remove_vals": []
    }
    costs = {}

    for key in distributions["properties"]:

        values = sorted([int(i) for i in distributions[key]["distribution"].keys()])
        rang = distributions[key]["max_val"] - distributions[key]["pruning_step"] * distributions[key]["pruning_depth"]
        num_removed = 0
        removed_vals = []

        if rang < distributions[key]["min_val"]:
            continue

        while values[-1] > rang:
            remove_val = str(values.pop())
            num_removed += distributions[key]["distribution"][remove_val]
            removed_vals.append(remove_val)

        costs[key] = {
            "num_removed": num_removed,
            "remove_vals": removed_vals
        }

        if (num_removed >= max_cost["num_removed"]) and (total_removed + num_removed <= num_to_remove):
            max_cost = {
                "key": key,
                "num_removed": num_removed,
                "remove_vals": removed_vals
            }

    return max_cost


def _greedy_balanced_removal(distributions, num_candidates, temp_dir):
    num_to_remove = num_candidates - TASK["candidates"]["prepare"]["max_molecules"]
    logger.info("Targeting to remove %d from %d candidates",
                num_candidates, num_to_remove)
    total_removed = 0
    step = 0
    max_cost_was_none = False
    while True:
        max_cost = _select_property_values_for_removal(distributions, total_removed, num_to_remove)

        if max_cost["key"] is None:
            # there is no full category can be included without surpassing the num_to_remove
            if max_cost_was_none:
                logger.debug("Cannot remove further full step from any property and further updating not helpful."
                             " Removing %d candidates.", num_to_remove - total_removed)
                return _remove_remaining_with_largest_values(distributions, num_to_remove - total_removed, temp_dir)
            else:
                logger.debug("Cannot remove further full step from any property, updating the source file.")
                remaining_candidates = _update_distributions_source_file(distributions)
                total_removed = num_candidates - remaining_candidates
                logger.info("%d candidates remain >> true number of removed = %d", remaining_candidates, total_removed)
                max_cost_was_none = True
                continue
        else:
            max_cost_was_none = False

        if max_cost["num_removed"] == 0:
            logger.debug("Nothing to remove >> increasing the pruning depth")
            for key in distributions["properties"]:
                distributions[key]["pruning_depth"] += 1
            continue

        step += 1

        if len(max_cost["remove_vals"]) > 0:
            for key in distributions["properties"]:
                if max_cost["key"] == key:
                    logger.debug("Removing values from %s distributions", key)
                    for i in max_cost["remove_vals"]:
                        del distributions[key]["distribution"][i]
                else:
                    distributions[key]["pruning_depth"] += 1

        total_removed += max_cost["num_removed"]
        logger.info("In %d step of pruning %d candidates were selected for removal based on property %s. "
                    "Keeping in total %d candidates", step, max_cost["num_removed"], max_cost["key"],
                    num_candidates - total_removed)


def _download_files_with_values(files, directory):
    ftp = Ftp("ftp.ncbi.nlm.nih.gov")
    for file in files:
        ftp.mirror_directory("pubchem/RDF/descriptor/compound",
                             ".*" + file + ".*",
                             directory)


def _create_filter_rdf_tasks(filters, in_dir, out_dir):
    output = []
    for file in os.listdir(in_dir):
        for key, filter_ in filters.items():
            if key not in file:
                continue
            output_file_name = file.replace(".ttl.gz", ".txt.gz")
            output.append({
                "property": key,
                "filter": filter_,
                "input": os.path.join(in_dir, file),
                "output": os.path.join(out_dir, output_file_name)
            })
    return output


def _filter_rdf_files(task):
    """Values lines are in format '.*has-value.*'."""
    stats = collections.defaultdict(int)
    filter_ = _filter_function_factory(task["filter"], stats)
    if os.path.exists(task["output"]):
        logger.info("Skipping: %s", task["input"])
        return
    logger.info("Filtering: %s", task["input"])
    with WritableValueFile(task["output"], stats) as value_file:
        with PubchemTurtleFile(task["input"]) as turtle:
            for cid, value in turtle:
                if filter_(value):
                    value_file.consume(cid, value)
    logger.info("Finished filtering: %s", task["input"])


def _filter_function_factory(filter_, stats):
    def range_filter(value):
        return filter_["value"][0] <= float(value) <= filter_["value"][1]

    if filter_["type"] == "range":
        return range_filter

    if filter_["type"] == "structure":
        return _create_structure_filter_function(filter_, stats)


def _create_structure_filter_function(filter_, stats):
    smiles_filter, mol_filter, smiles_atoms = _create_substructure_filter(
        Chem.MolFromSmarts(filter_["smiles"]), use_chirality=True)

    logger.debug("Used molecule SMILES:")
    logger.debug("   %s", filter_["smiles"])

    logger.debug("Used Fragments SMILES:")
    frag_smiles_filters = []
    frag_mol_filters = []
    for mol in _prepare_fragments_molecules(filter_):
        logger.debug("   %s", Chem.MolToSmiles(mol))
        frag_smiles_filter, frag_mol_filter, _ = _create_substructure_filter(
            mol, ignored_atoms=smiles_atoms)
        frag_smiles_filters.append(frag_smiles_filter)
        frag_mol_filters.append(frag_mol_filter)

    def substructure_filter(smiles):

        if not smiles_filter(smiles):
            stats["smiles_filter"] += 1
            return False

        frag_smiles_res = [pred(smiles) for pred in frag_smiles_filters]
        if not any(frag_smiles_res):
            stats["frag_smiles_filter"] += 1
            return False

        molecule = Chem.MolFromSmiles(smiles, sanitize=False)
        if molecule is None:
            logger.error("Invalid SMILES: %s", smiles)
            return False

        try:
            # As we do not use sanitize in loading we need to do some
            # here. It is also required for AddHs.
            molecule.UpdatePropertyCache()
        except ValueError as ex:
            stats["invalid"] += 1
            logger.error("Invalid molecule: %s : %s", smiles, ex)
            return False

        # It is actually faster with chirality.
        if not mol_filter(molecule):
            stats["mol_filter"] += 1
            return False

        # We add Hydrogens here as the SMILES filter is not using them
        # but the fragment filter can.
        molecule = Chem.AddHs(molecule)

        if not any([should_test and pred(molecule)
                    for pred, should_test in
                    zip(frag_mol_filters, frag_smiles_res)]):
            stats["frag_mol_filter"] += 1
            return False

        atoms = set([atom.GetSymbol() for atom in molecule.GetAtoms()])
        for atom in atoms:
            if atom not in CHEM["allowed_atoms"]:
                stats["prohibited_atom_filter"] += 1
                return False

        stats["pass"] += 1
        return True

    return substructure_filter


def _create_substructure_filter(sub_mol, use_chirality=False, ignored_atoms=None):
    if ignored_atoms is None:
        ignored_atoms = []

    # 0 is for * 1 for H
    required_atoms = set([
        atom.GetSymbol().lower()
        for atom in sub_mol.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom not in ignored_atoms
    ])

    def smiles_filter(smiles):
        # We use lower case, this will give us some false positives
        # but should prevent any false negatives.
        smiles = smiles.lower()
        for atom in required_atoms:
            if atom not in smiles:
                return False
        return True

    def molecule_filter(molecule):
        return molecule.HasSubstructMatch(sub_mol, useChirality=use_chirality)

    return smiles_filter, molecule_filter, required_atoms


def _prepare_fragments_molecules(filter_):
    expanded_fragment = expand_fragment(
        filter_["fragment"],
        [
            os.path.join(filter_["substrates"], file)
            for file in os.listdir(filter_["substrates"])
            if file.lower().endswith(".sdf")
        ],
        ignore_atom_types=True)
    fragments_molecules = fragment_to_rdkit_molecules(expanded_fragment)
    # When we create molecules above, they have the right SMILES but
    # invalid SMARST for example we get:
    # SMILES: [*]C([*])([*])Br
    # SMARTS: *-,:C(-,:*)(-,:*)-,:[Br]
    # However we SMARTS should be: [#35]-[#6](-[*])(-[*])-[*]
    # A workaround seems to be convert the molecule
    # to SMILES and then use the string as a SMARTS.
    # TODO Find cause and remove workaround.
    fragments_molecules = [
        Chem.MolFromSmarts(Chem.MolToSmiles(mol))
        for mol in fragments_molecules
    ]
    return fragments_molecules


def _convert_rdkit_molecule_to_molecule(rdkit_mol):
    class Molecule(object):

        def __init__(self, atoms, bonds):
            self._atoms = atoms
            self._bonds = bonds

        def atoms(self):
            return self._atoms

        def get_bond(self, left_atom, right_atom):
            origin = left_atom.idx
            target = right_atom.idx

            if origin > target:
                origin, target = target, origin

            for bond in self._bonds:
                if bond.origin == origin and bond.target == target:
                    return bond
            return None

    class Atom(object):

        def __init__(self, rdkit_atom):
            self.idx = rdkit_atom.GetIdx()
            self.atom_type = rdkit_atom.GetSymbol()

    class Bond(object):

        def __init__(self, rdkit_bond):
            self.origin = rdkit_bond.GetBeginAtomIdx()
            self.target = rdkit_bond.GetEndAtomIdx()
            if self.origin > self.target:
                self.origin, self.target = self.target, self.origin
            self.type = rdkit_bond.GetBondType()

    atoms = [Atom(atom) for atom in rdkit_mol.GetAtoms()]
    atoms.sort(key=lambda atom: atom.idx)
    bonds = [Bond(bond) for bond in rdkit_mol.GetBonds()]
    return Molecule(atoms, bonds)


def _create_tasks(temp_dir, filter_tasks, step_size):
    files_info = []
    properties = set()
    for task in filter_tasks:
        properties.add(task["property"])
        task_info = read_json(task["output"] + ".json")
        task_info["property"] = task["property"]
        task_info["file"] = task["output"]
        files_info.append(task_info)
    min_cid = min([info["min"] for info in files_info])
    max_cid = max([info["max"] for info in files_info])
    logger.info("CID Range: %d, %d", min_cid, max_cid)

    max_cid_len = len(str(max_cid))
    output = []
    for interval_start in range(min_cid, max_cid, step_size):
        interval_end = min(interval_start + step_size - 1, max_cid)
        new_task = {
            "start": interval_start,
            "end": interval_end,
            "files": {},
            "output": os.path.join(
                temp_dir,
                "filter-{}-{}.txt.gz".format(
                    str(interval_start).zfill(max_cid_len),
                    str(interval_end).zfill(max_cid_len)))
        }
        for property_ in properties:
            files_for_property = []
            for info in files_info:
                if info["property"] != property_:
                    continue
                if info["min"] > interval_end:
                    continue
                if info["max"] < interval_start:
                    continue
                files_for_property.append(info["file"])
            new_task["files"][property_] = files_for_property
        output.append(new_task)
    return output


def _merge_value_files(task):
    if os.path.exists(task["output"]):
        logger.info("Skipping: %s", task["output"])
        return
    logger.info("Preparing merging for: %s", task["output"])
    properties = {}
    for property_, files in task["files"].items():
        values = []
        for file in files:
            logger.debug("Reading: %s", file)
            with ReadableValueFile(file) as reader:
                # We can store values here as well, but that would cost
                # lot of extra memory (especially for SMILES).
                values.extend([key for key, value in reader])
        properties[property_] = sorted(values)

    logger.debug("Merging data ...")
    counter = 0
    with gzip.open(task["output"], "wt") as stream:
        for value in _iterate_common_values(
                properties.values(), task["start"], task["end"]):
            stream.write("{}\n".format(value))
            counter += 1
    with open(task["output"] + ".json", "w") as stream:
        json.dump({"count": counter}, stream)


def _iterate_common_values(lists, start, end):
    sources = []
    for item in lists:
        generator = itertools.chain(item)
        value = next(generator)
        # Skip all values until start.
        while value < start:
            value = next(generator)
        sources.append({"generator": generator, "value": value})

    try:
        while True:
            values = [source["value"] for source in sources]
            min_value = min(values)
            max_value = max(values)
            if max_value > end:
                # One of the values reached the end.
                return
            if min_value == max_value:
                yield min_value
            for source in sources:
                if source["value"] == min_value:
                    source["value"] = next(source["generator"])
    except StopIteration:
        pass


# TODO Move to extra file.
BOND_TYPE_TO_ORDER = {
    1: Chem.rdchem.BondType.SINGLE,
    2: Chem.rdchem.BondType.DOUBLE
}


def fragment_to_rdkit_molecules(fragment):
    return [
        _fragment_to_rdkit_molecule(fragment)
        for fragment in _unpack_fragment(fragment)
    ]


def _fragment_to_rdkit_molecule(fragment):
    # TODO Add types from SDF file.
    molecule = Chem.RWMol(Chem.Mol())
    for type_ in fragment["types"]:
        molecule.AddAtom(Chem.Atom(type_))
    for bond in fragment["bonds"]:
        molecule.AddBond(bond[0], bond[1], BOND_TYPE_TO_ORDER[bond[2]])
    return molecule


def _unpack_fragment(fragment):
    """Expand atom symbol groups."""
    fragments = [fragment]
    for index in range(0, len(fragment["types"])):
        types = _unpack_atom_symbol(fragment["types"][index])
        fragments = [
            {
                "names": f["names"],
                "types": f["types"][:index] + [type_] + f["types"][index + 1:],
                "bonds": f["bonds"]
            }
            for f in fragments
            for type_ in types
        ]
    return fragments


def _unpack_atom_symbol(symbol):
    if symbol == "X":
        return ["F", "Cl", "Br", "I", "At"]
    else:
        return [symbol]


# region Expand fragment

FragmentExpansion = collections.namedtuple(
    "Candidate", ["atom_type", "sources", "bonds"])

# TODO Move to extra file.
BOND_ORDER_TO_TYPE = {
    Chem.rdchem.BondType.SINGLE: 1,
    Chem.rdchem.BondType.DOUBLE: 2
}


def expand_fragment(fragment_to_expand, files, ignore_atom_types=False):
    molecules = _load_sdf_files(files)
    fragments = fragment_to_rdkit_molecules(fragment_to_expand)
    molecules_matches = []
    for molecule in molecules:
        matches = []
        for fragment in fragments:
            matches.extend(molecule.GetSubstructMatches(fragment))
        molecules_matches.append(matches)

    # For each molecule store array of expansions.
    # Expansions are defined for each math - ie. the can not be mixed together
    # as they corresponds to different part of the molecule.
    expansions = []
    for index, matches in enumerate(molecules_matches):
        # Use set so we keep only unique items.
        molecules_expansions = {
            tuple(_propose_expansions(molecules[index], match))
            for match in matches
        }
        expansions.append(molecules_expansions)

    # Optionally we can remove atom types.
    if ignore_atom_types:
        expansions = _replace_atom_type_with_wildcard(expansions)

    # Pick such tuples that are in every molecule.
    output = []
    candidates = set([item for item in itertools.chain(*expansions)])
    for candidate in candidates:
        if not _id_supported_by_all_molecules(candidate, expansions):
            continue
        output.append(
            _expand_fragment_from_candidate(fragment_to_expand, candidate))
    assert len(output) == 1, "Only one expanded fragment is supported."
    return output[0]


def _load_sdf_files(files):
    output = []
    for file in files:
        suppl = Chem.SDMolSupplier(file)
        output.extend([Chem.AddHs(mol) for mol in suppl if mol is not None])
    return output


def _propose_expansions(molecule, fragment):
    """
    Return a list of expansion candidates.

    Each candidate is identified by atom type, and all connections to fragment.
    """
    fragment_atoms = set(fragment)

    new_atoms = {}
    for bond in molecule.GetBonds():
        bond_atoms = {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
        intersection = fragment_atoms & bond_atoms
        intersection_size = len(intersection)
        if intersection_size == 2:
            # Bond is already in the fragment
            continue
        elif intersection_size == 1:
            # Bond can expand the fragment
            begin_index = fragment.index(intersection.pop())
            end_idx = (bond_atoms - fragment_atoms).pop()
            bond_type = bond.GetBondType()
            if end_idx not in new_atoms:
                new_atoms[end_idx] = {
                    "type": molecule.GetAtomWithIdx(end_idx).GetSymbol(),
                    "sources": [],
                    "bonds": []
                }
            new_atoms[end_idx]["sources"].append(begin_index)
            new_atoms[end_idx]["bonds"].append(bond_type)

        else:
            # Bond is outside the fragment.
            pass

    return sorted([
        FragmentExpansion(item["type"], tuple(item["sources"]),
                          tuple(item["bonds"]))
        for item in new_atoms.values()
    ])


def _replace_atom_type_with_wildcard(expansions):
    return [[
        tuple([FragmentExpansion("*", *item[1:]) for item in expansion])
        for expansion in mol
    ] for mol in expansions]


def _id_supported_by_all_molecules(expansion, molecules_expansions):
    """Return True if expansion is supported by all molecules."""

    def is_subset(subset, superset):
        for key, count in subset.items():
            if key not in superset:
                return False
            if count > superset[key]:
                return False
        return True

    for molecule in molecules_expansions:
        if expansion in molecule:
            continue

        # Check for subset.
        expansion_counter = collections.Counter(expansion)
        if any([is_subset(expansion_counter, collections.Counter(item))
                for item in molecule]):
            continue

        return False
    return True


def _expand_fragment_from_candidate(fragment, candidate):
    new_types = \
        fragment["types"] + \
        [candidate.atom_type for candidate in candidate]
    new_names = \
        [type_ + str(index + 1) for index, type_ in enumerate(new_types)]

    # Use this to create new instance of array.
    new_bonds = [] + fragment["bonds"]
    for index, item in enumerate(candidate):
        # Move by offset.
        index += len(fragment["types"])
        for source, bond in zip(item.sources, item.bonds):
            new_bonds.append((index, source, BOND_ORDER_TO_TYPE[bond]))

    return {
        "names": new_names,
        "types": new_types,
        "bonds": new_bonds
    }


# endregion

# region Download from PubChem by CID

def _count_candidates(temp_dir):
    cids = {}
    counter = 0
    for filter_file in os.listdir(temp_dir):
        if "filter" not in filter_file or not filter_file.endswith(".txt.gz") or "subdir" in filter_file:
            continue

        filter_path = os.path.join(temp_dir, filter_file)
        with gzip.open(filter_path, "rt") as in_stream:
            for line in in_stream:
                cid = line.rstrip()
                subset_id = str(int(cid) // SUBSET_SIZE)
                if subset_id not in cids.keys():
                    cids[subset_id] = []
                cids[subset_id].append(cid)
                counter += 1

    return counter, cids


def _update_input_files_for_download(temp_dir, output_dir, candidates_to_keep):
    logger.debug("Mapping local directory structure.")
    cids = []
    total_molecules = 0

    for filter_file in sorted(os.listdir(temp_dir)):
        if "filter" not in filter_file or not filter_file.endswith(".txt.gz") or "subdir" in filter_file:
            continue
        filter_path = os.path.join(temp_dir, filter_file)
        with gzip.open(filter_path, "rt") as in_stream:
            for line in in_stream:
                cid = line.rstrip()
                subset_id = str(int(cid) // SUBSET_SIZE)
                if subset_id in candidates_to_keep.keys() and cid in candidates_to_keep[subset_id]:
                    cids.append(int(cid))

    batch_size = math.ceil(len(cids) / ENV["num_parallel_cpu"])
    if batch_size > ENV["pubchem"]["batch_size"]:
        batch_size = ENV["pubchem"]["batch_size"]

    num_batches = math.ceil(len(cids) / batch_size)
    logger.debug("Assigning candidates to %d batches of size %d for download", num_batches, batch_size)
    for num in range(0, num_batches):
        logger.debug("Processing %d batch", num + 1)
        out_file = os.path.join(temp_dir, "filter-{}_subdir.txt.gz".format(num))

        if os.path.exists(out_file + ".json"):
            info = read_json(out_file + ".json")
            total_molecules += info["count"]
            continue

        with gzip.open(out_file, "wt") as out_stream:
            molecule_counter = 0
            while (molecule_counter <= batch_size) and len(cids) > 0:
                cid = cids.pop(0)
                molecule_counter += 1
                dir_name = "subdir{:=04d}".format(math.ceil(
                    float(total_molecules + molecule_counter) / ENV["max_candidates_in_dir"]))
                out_stream.write("{},{}\n".format(cid, dir_name))

            total_molecules += molecule_counter

        with open(out_file + ".json", "w") as out_json_stream:
            json.dump({"count": molecule_counter}, out_json_stream)

    num_subdirs = math.ceil(float(total_molecules) / ENV["max_candidates_in_dir"])
    logger.debug("Creating %d subdirs for candidates in %s", num_subdirs, output_dir)
    for i in range(1, num_subdirs + 1):
        sub_dir = os.path.join(output_dir, "subdir{:=04d}".format(i))
        os.makedirs(sub_dir, exist_ok=True)

    return total_molecules


def _create_download_tasks(temp_dir, output_dir):
    download_tasks = []

    for filter_file in os.listdir(temp_dir):
        if "filter" not in filter_file or not filter_file.endswith("_subdir.txt.gz"):
            continue
        download_tasks.append({
            # Pick one of the query windows.
            "input": os.path.join(temp_dir, filter_file),
            "output": os.path.join(output_dir, "sdf"),
            "temp_dir": os.path.join(paths.temp(), "pubchem_download")
        })

    return download_tasks


class PubchemSdfBatchDownloader(object):
    """
    Download SDF with batch of molecules with given CDF from PubChem.
    """

    def __init__(self, can_query=None):
        if can_query is None:
            can_query = lambda: True
        self._last_request_xml = None
        self._can_query = can_query

    def download(self, batch, output_path):
        logger.debug("Requesting data from PubChem ...")

        response = self._http_post(self._request_xml(batch))
        if not self._are_data_ready(response):
            response = self._wait_till_data_are_ready(response)

        download_url = self._get_download_url(response.text)
        self._download_from_ftp(download_url, output_path)

    @staticmethod
    def _request_xml(cids):
        request_xml_molecules = "\n".join([
            "<PCT-ID-List_uids_E>" + cid + "</PCT-ID-List_uids_E>"
            for cid in cids
        ])
        return """
<?xml version="1.0"?>
<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">
<PCT-Data>
  <PCT-Data_input>
    <PCT-InputData>
      <PCT-InputData_download>
        <PCT-Download>
          <PCT-Download_uids>
            <PCT-QueryUids>
              <PCT-QueryUids_ids>
                <PCT-ID-List>
                  <PCT-ID-List_db>pccompound</PCT-ID-List_db>
                  <PCT-ID-List_uids>
                  """ + request_xml_molecules + """
                  </PCT-ID-List_uids>
                </PCT-ID-List>
              </PCT-QueryUids_ids>
            </PCT-QueryUids>
          </PCT-Download_uids>
          <PCT-Download_format value="sdf"/>
          <PCT-Download_compression value="gzip"/>
          <PCT-Download_use-3d value="true"/>
        </PCT-Download>
      </PCT-InputData_download>
    </PCT-InputData>
  </PCT-Data_input>
</PCT-Data>
"""

    def _http_post(self, xml):

        # Wait until we can download the data.
        while not self._can_query():
            time.sleep(0.6)

        self._last_request_xml = xml
        headers = {"Content-Type": "application/xml"}
        response = requests.post(
            "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi",
            data=xml, headers=headers, timeout=10 * 60)
        if not response.status_code == 200:
            raise RuntimeError(
                "Request failed.\n    Code: {}\n    Response: {}".format(
                    response.status_code, response.text))
        return response

    @staticmethod
    def _are_data_ready(response):
        return "<PCT-Status value=\"success\"/>" in response.text

    def _wait_till_data_are_ready(self, response):
        body = response.text
        start_tag = "<PCT-Waiting_reqid>"
        end_tag = "</PCT-Waiting_reqid>"
        try:
            reqid = body[
                    body.index(start_tag) + len(start_tag):
                    body.index(end_tag)
                    ]
        except ValueError:
            raise RuntimeError(
                "Unexpected response: {}\n{}\n\nLast request:\{}".format(
                    response.status_code, body, self._last_request_xml))

        while "<PCT-Status value=\"running\"/>" in response.text:
            logger.debug("Waiting for PubChem to prepare data ...")
            time.sleep(30)
            response = self._http_post(self._wait_xml(reqid))

        if not self._are_data_ready(response):
            raise RuntimeError(
                "Request finished unsuccessfully.\n    Code: {}\n    Response: {}".format(
                    response.status_code, response.text))

        logger.debug("Data are ready")
        return response

    @staticmethod
    def _wait_xml(reqid):
        return """
<?xml version="1.0"?>
<!DOCTYPE PCT-Data PUBLIC "-//NCBI//NCBI PCTools/EN" "http://pubchem.ncbi.nlm.nih.gov/pug/pug.dtd">
<PCT-Data>
<PCT-Data_input>
<PCT-InputData>
  <PCT-InputData_request>
    <PCT-Request>
      <PCT-Request_reqid>""" + reqid + """</PCT-Request_reqid>
         <PCT-Request_type value="status"/>
       </PCT-Request>
     </PCT-InputData_request>
   </PCT-InputData>
 </PCT-Data_input>
</PCT-Data>
"""

    def _get_download_url(self, body):
        start_tag = "<PCT-Download-URL_url>"
        end_tag = "</PCT-Download-URL_url>"
        return body[
               body.index(start_tag) + len(start_tag):
               body.index(end_tag)
               ]

    def _download_from_ftp(self, url, output_path):
        if not url.startswith("ftp://"):
            raise RuntimeError(
                "Invalid protocol in '" + url + "' FTP expected.")

        download_url = url[len("ftp://"):]
        domain = download_url[:download_url.index("/")]
        path = download_url[download_url.index("/") + 1:]

        ftp = Ftp(domain)
        ftp.mirror_file(path, output_path)


def _download_from_pubchem_by_cid(task):
    with open(task["input"] + ".json", "r") as stream:
        info = json.load(stream)

    logger.info("Downloading molecules for: %s", task["input"])

    counter = 0
    counter_step = 100000
    counter_next = counter_step

    downloader = _create_pubchem_sdf_downloader()

    batch_name = os.path.basename(task["input"])[:-len(".txt.gz")]

    os.makedirs(task["temp_dir"], exist_ok=True)
    for batch, name in _iterate_cid_of_missing_in_batches(task):
        temp_file = os.path.join(task["temp_dir"], batch_name + "_" + name + ".sdf.gz")

        if len(batch) > 0:
            if not os.path.exists(temp_file):
                downloader.download(batch.keys(), temp_file)

            _split_sdf_file(temp_file, task, batch)

            counter += len(batch)
            if counter > counter_next:
                logger.info("Downloading batch %s %d/%d",
                            batch_name, counter, info["count"])
                counter_next += counter_step


def _assign_query_time(entrez=False, proc_per_second=2):
    # Define time windows in which each thread can query a database.
    # This is just a simple solution but it should be ok.

    process = multiprocessing.current_process()
    process_id = process._identity[0]
    window_count = min(ENV["num_parallel_cpu"], 60)
    window_size = math.floor(60 / window_count)

    if entrez:
        if proc_per_second > 3:
            raise RuntimeError("No more than 3 requests per a second are allowed for Entrez querying.")
        query_windows = []
        for x in range(0, window_count // proc_per_second):
            query_windows.append(
                [int((x + window_count // proc_per_second * y)) for y in range(0, window_size * proc_per_second)])
        # Pick query window used in this thread.
        query_window = query_windows[process_id % len(query_windows)]
        logger.debug("Using query times from: " + str(query_window))
        return query_window
    else:
        query_windows = [(x, x + window_size - 1) for x in range(0, 60, window_size)]
        # Pick query window used in this thread.
        query_window = query_windows[process_id % len(query_windows)]
        window_start = query_window[0]
        window_end = query_window[1]
        logger.debug("Using query window: %d - %d", window_start, window_end)
        return window_start, window_end


def _create_pubchem_sdf_downloader():
    window_start, window_end = _assign_query_time()

    def is_in_query_window():
        return window_start <= datetime.now().second <= window_end

    return PubchemSdfBatchDownloader(is_in_query_window)


def _iterate_cid_of_missing_in_batches(task):
    batch_size = ENV["pubchem"]["batch_size"]
    counter = 0
    with gzip.open(task["input"], "rt") as stream:
        molecules = {}
        for line in stream:
            counter += 1
            cid, subdir = line.rstrip().split(",")
            sdf_path = os.path.join(task["output"], subdir, cid + ".sdf")
            if not os.path.exists(sdf_path):
                molecules[cid] = subdir
            if len(molecules.keys()) > batch_size:
                yield molecules, str(counter)
                molecules = {}
        if len(molecules.keys()) > 0:
            yield molecules, str(counter)


def _split_sdf_file(input_path, task, batch):
    logger.debug("Extracting SDF files from %s", input_path)

    def _process_failed_file(file_path, err_msg):
        if os.path.exists(file_path):
            logger.error("Following error occurred %s, when processing downloaded file %s", err_msg, file_path)
            failed_dir = os.path.join(paths.temp(), "pubchem_download", "failed")
            os.makedirs(failed_dir, exist_ok=True)
            failed_path = os.path.join(failed_dir, os.path.basename(file_path)) + "_" + str(time.time())
            os.rename(file_path, failed_path)

    output_stream = None
    try:
        with gzip.open(input_path, "rt") as input_stream:
            for line in input_stream:
                if line.startswith("$$$$"):
                    output_stream.close()
                    output_stream = None
                    continue
                if output_stream is not None:
                    output_stream.write(line)
                    continue
                if line.isspace():
                    continue
                # New molecule and line contains their name.
                name = line.rstrip()
                output_path = os.path.join(task["output"], batch[name], name + ".sdf")
                output_stream = open(output_path, "w")
                output_stream.write(line)
    except OSError as os_err:
        _process_failed_file(input_path, os_err)
    except EOFError as eof_err:
        _process_failed_file(input_path, eof_err)
    finally:
        if output_stream is not None:
            output_stream.close()
        logger.debug("Extracting SDF files from %s done", input_path)
        remove_existing(input_path)


def _evaluate_sdf_output_dir(output_dir):
    molecule_counter = 0
    for subdir in os.listdir(output_dir):
        if "subdir" not in subdir:
            continue
        subdir_path = os.path.join(output_dir, subdir)
        molecule_counter += len(os.listdir(subdir_path))

    return molecule_counter


# endregion


def workaround_to_investigate():
    # TODO Investigate difference in SMARTS as the SMILES is the same.

    mol_smiles = "C1CNCCC1C2=NC3=CC(=C(C=C3N2CC4=CC=C(C=C4)CBr)Cl)Cl"
    mol = Chem.AddHs(Chem.MolFromSmiles(mol_smiles))

    frag_smiles = "[*]C([*])([*])Br"
    frag = Chem.MolFromSmarts(frag_smiles)
    # frag.UpdatePropertyCache()

    print(Chem.MolToSmiles(mol))
    print(Chem.MolToSmarts(frag))
    print(mol.HasSubstructMatch(frag))

    print("-" * 80)
    for atom in frag.GetAtoms():
        print("ATOM", atom.GetIdx(), atom.GetSymbol())
    for bond in frag.GetBonds():
        print("BOND",
              bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())

    print("\n" * 6)

    mol_smiles = "C1CNCCC1C2=NC3=CC(=C(C=C3N2CC4=CC=C(C=C4)CBr)Cl)Cl"
    fragment = {
        'types': ['Br', 'C', '*', '*', '*'],
        'names': ['X1', 'C2', '*3', '*4', '*5'],
        'bonds': [(1, 0, 1), (2, 1, 1), (3, 1, 1), (4, 1, 1)]
    }
    fragments_molecules = fragment_to_rdkit_molecules(fragment)

    mol = Chem.AddHs(Chem.MolFromSmiles(mol_smiles))
    print(Chem.MolToSmiles(mol))

    for frag_mol in fragments_molecules:
        frag_mol.UpdatePropertyCache()

        print(Chem.MolToSmarts(frag_mol))
        print(Chem.MolToSmiles(frag_mol))
        # Workaround to fix SMARTS.
        frag_mol_smiles = Chem.MolFromSmarts(Chem.MolToSmiles(frag_mol))

        print(mol.HasSubstructMatch(frag_mol))
        print(mol.HasSubstructMatch(frag_mol_smiles))

        print("-" * 80)
        for atom in frag_mol.GetAtoms():
            print("ATOM", atom.GetIdx(), atom.GetSymbol())
        for bond in frag_mol.GetBonds():
            print("BOND",
                  bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
                  bond.GetBondType())
