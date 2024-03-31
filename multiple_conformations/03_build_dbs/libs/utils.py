#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import re
import os
import gzip
import json
import urllib.request
import shutil
import logging
import datetime
from libs.configuration import get_num_cpu, ENV
from libs import paths
from libs.database import pdb as pdb_db
from external import clustal_w as clustal_w_service


logger = logging.getLogger(__name__)


class ObjectCache(object):
    """
    Cache objects into files.

    The manager must implement following methods:
    * create(self, id)
    * to_json(self, object)
    * from_json(self, json)
    and must has static property:
    * cache_name
    """

    def __init__(self, manager, identifier):
        self.id = identifier
        self.manager = manager
        self.stream = None
        self.cache_directory = os.path.join(ENV["cache_directory"], manager.cache_name)
        os.makedirs(self.cache_directory, exist_ok=True)

    def __enter__(self):
        if os.path.exists(self.path()):
            return self.read_from_cache()
        else:
            return self.cache_object()

    def __exit__(self, *args):
        pass

    def path(self):
        return os.path.join(self.cache_directory, self.id)

    def cache_object(self):
        path = self.path()
        object_value = self.manager.create(self.id)
        with open(path, "w") as stream:
            json.dump(self.manager.to_json(object_value), stream)
        return object_value

    def read_from_cache(self):
        path = self.path()
        with open(path, "r") as stream:
            value = json.load(stream)
        return self.manager.from_json(value)


class UrlCache(object):
    """
    Cache HTTP GET request.

    TODO Add remove need to specify name and enable use of arbitrary URL.
    """

    def __init__(self, identifier, url, cache_name="url-cache", mode="r"):
        self.identifier = str(identifier)
        self.url = url
        self.cache_directory = os.path.join(ENV["cache_directory"], cache_name)
        os.makedirs(self.cache_directory, exist_ok=True)
        self.mode = mode
        self.stream = None

    def __enter__(self):
        if not os.path.exists(self._path()):
            self.download()
        self.stream = open(self._path(), self.mode, encoding="UTF-8")
        return self.stream

    def __exit__(self, *args):
        self.stream.close()

    def _path(self):
        return os.path.join(self.cache_directory, self.identifier)

    def as_path(self):
        if not os.path.exists(self._path()):
            self.download()
        return os.path.join(self.cache_directory, self.identifier)

    def download(self):
        with urllib.request.urlopen(self.url) as response:
            with open(self._path(), "wb") as output_stream:
                output_stream.write(response.read())


def dict_append(dictionary, key, value):
    list_to_add_to = dictionary.get(key, None)
    if list_to_add_to is None:
        list_to_add_to = []
        dictionary[key] = list_to_add_to
    list_to_add_to.append(value)


def fetch_url(identifier, url, cache_name):
    cache_directory = os.path.join(ENV["cache_directory"], cache_name)
    os.makedirs(cache_directory, exist_ok=True)
    local_path = os.path.join(cache_directory, str(identifier))
    if os.path.exists(local_path):
        return local_path
    else:
        try:
            with urllib.request.urlopen(url) as response:
                with open(local_path, "wb") as output_stream:
                    output_stream.write(response.read())
        except urllib.error.HTTPError:
            raise FileNotFoundError
    return local_path


def copy_files(paths, output_directory):
    for path in paths:
        file_name = os.path.basename(path)
        shutil.copy(path, os.path.join(output_directory, file_name))


def remove_existing(path):
    if os.path.exists(path):
        os.remove(path)


def blast_uniprot(url, output_path):
    with urllib.request.urlopen(url) as response:
        with open(output_path, "wb") as output_stream:
            output_stream.write(response.read())
    return os.path.exists(output_path)


def _compute_parallel_set_sizes(batch_size, set_size=None):
    num_processess = get_num_cpu(ENV["num_parallel_cpu"])
    set_sizes = [batch_size // num_processess] * num_processess

    if batch_size % num_processess == 0:
        return set_sizes
    else:
        for i in range(0, batch_size % num_processess):
            set_sizes[i] += 1
        return set_sizes


def compute_range_for_sets(batch_size, compute_size_fnc=None, set_size=None):
    begin_ends = []
    begin = 0
    if compute_size_fnc is None:
        compute_size_fnc = _compute_parallel_set_sizes
    for size in compute_size_fnc(batch_size, set_size):
        end = begin + size
        begin_ends.append((begin, end))
        begin = end

    return begin_ends


def split_fragments_to_sets(fragments):
    fragments_sets = []

    for fragment_range in compute_range_for_sets(len(fragments)):
        fragments_sets.append(fragments[fragment_range[0]:fragment_range[1]])

    return fragments_sets


def split_files_in_subdirs_to_sets(in_dir):
    files = []
    for subdir_name in os.listdir(in_dir):
        subdir = os.path.join(in_dir, subdir_name)
        for file_name in os.listdir(subdir):
            files.append(os.path.join(subdir, file_name))

    files_sets = []
    for files_range in compute_range_for_sets(len(files)):
        files_sets.append(files[files_range[0]:files_range[1]])

    return files_sets


_CONSOLE_HANDLER = None


def init_logging():
    global _CONSOLE_HANDLER
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        fmt="%(asctime)s [%(levelname)s] %(thread)d - %(message)s",
        datefmt="%H:%M:%S")
    handler.setFormatter(formatter)
    handler.setLevel(logging.INFO)
    logging.getLogger().addHandler(handler)
    logging.getLogger().setLevel(logging.NOTSET)
    _CONSOLE_HANDLER = handler


def set_logging_level(level):
    _CONSOLE_HANDLER.setLevel(level)


def read_json(path):
    with open(path) as stream:
        return json.load(stream)


def read_task_info():
    with open(paths.task_info()) as in_stream:
        return json.load(in_stream)


def load_json_lines(path, gzipped=False):
    output = []
    if gzipped:
        with gzip.open(path) as in_stream:
            for line in in_stream:
                output.extend(json.loads(line.decode("utf-8")))
    else:
        with open(path) as in_stream:
            for line in in_stream:
                output.extend(json.loads(line))

    return output


def among_candidate_processing_statuses(values):
    for value in values:
        if is_candidate_processing_status(value):
            return True

    return False


def is_candidate_processing_status(value):
    """"
    possible states:
    initial                     - after _export_candidates_info() in 08_candidates_prepare.py
    sdf-downloaded              - after _query_pubchem_for_substrates() in 08_candidates_prepare.py
    pdbqt-prepared              - after _prepare_candidates_pdbqts() in 08_candidates_prepare.py
    docking-prepared            - after whole 08_candidates_prepare.py
    ready_for_transfer2remote   - after the local part of set_up_screening() in remote_docking.py
    transferred2remote          - after the remote part of set_up_screening() in remote_docking.py
    docking-completed           - after 09_candidates_docking.py
    docking-failed              - after 09_candidates_docking.py
    """
    status = None
    info_path = os.path.join(paths.candidates(), "info.json")
    if os.path.exists(info_path):
        status = read_json(info_path)["status"]

    if value is not None:
        return status == value
    else:
        if status is None:
            return True
        else:
            return False


def update_candidate_processing_status(status):
    info_path = os.path.join(paths.candidates(), "info.json")
    if os.path.exists(info_path):
        info = read_json(info_path)
        info["status"] = status
        with open(info_path, "w") as out_stream:
            json.dump(info, out_stream, indent=2)


def values_array_to_dict(array):
    dictionary = {}
    for item in array:
        dictionary[tuple(item["name"])] = item["val"]
    return dictionary


class FileLogging(object):
    def __init__(self, directory, module_path):
        module_name = os.path.basename(module_path)[:-3]

        self.old_stdout = sys.stdout
        self.file_stream = None
        self.stdout_path = os.path.join(directory, module_name + ".out")

        log_path = os.path.join(directory, module_name + ".log")
        self.handler = logging.FileHandler(log_path)
        self.handler.setLevel(logging.NOTSET)
        formatter = logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(name)s : %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S")
        self.handler.setFormatter(formatter)
        logging.getLogger().addHandler(self.handler)

    def __enter__(self):
        sys.stdout = self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_stream is not None:
            self.file_stream.close()

        sys.stdout = self.old_stdout
        logging.getLogger().handlers = [
            handler for handler in logging.getLogger().handlers
            if not handler == self.handler]
        self.handler.close()

    def write(self, message):
        self.old_stdout.write(message)
        if self.file_stream is None:
            self.file_stream = open(self.stdout_path, "w")
        self.file_stream.write(message)

    def flush(self):
        if self.file_stream is not None:
            self.file_stream.flush()


def transform_files(directory, fnc, input_pattern, backup_template):
    names_to_process = {}
    for file_name in sorted(os.listdir(directory)):
        name = file_name[:file_name.find(".json")]
        if not re.match(input_pattern, name):
            logger.debug("Skipping file: %s named %s", file_name, name)
            continue
        if name not in names_to_process.keys():
            names_to_process[name] = []

        names_to_process[name].append(os.path.join(directory, file_name))

        if backup_template is not None:
            shutil.copy(
                os.path.join(directory, file_name),
                os.path.join(directory, backup_template.format(name))
            )

    for name in names_to_process.keys():
        fnc(names_to_process[name], names_to_process[name].copy())


def calculate_similarity(query_pdbid, candidate_pdbid, chains):
    path = os.path.join(paths.temp(), "similarity-input.tmp")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if len(chains) > 1:
        logger.error("There is more than one chain involved in the active site. "
                     "This functionality is not supported yet.")
        raise RuntimeError
    with open(path, "w") as output_stream:
        query_fasta = _get_fasta_chain(query_pdbid, chains[0])
        candidate_fasta = _get_fasta_chain(candidate_pdbid, chains[0])
        output_stream.write(query_fasta)
        output_stream.write(candidate_fasta)
    clustal_w_service.execute(path)
    path = os.path.join(paths.temp(), "similarity.log")
    chunks = [-1.0]
    with open(path, "r") as stream:
        for line in stream:
            if "Sequences (1:2) Aligned. Score:" in line:
                chunks = line.strip().split()
                break
    logger.debug("Removing similarity logfile: %s", path)
    os.remove(path)
    return float(chunks[-1])


def _get_fasta_chain(pdbid, chain):
    original_fasta = pdb_db.fetch_fasta(pdbid)
    pdb_fasta = ""
    with open(original_fasta, "r") as stream:
        line_s = stream.readline()
        while line_s:
            if line_s.startswith(">"):
                if _get_chain(line_s.strip()) == chain:
                    pdb_fasta += ">%4s:%1s|%s\n" %(pdbid, chain, str(datetime.datetime.now().timestamp()))
                    line_s = stream.readline()
                    while (not line_s.startswith(">")) and line_s:
                        pdb_fasta += line_s
                        line_s = stream.readline()
            line_s = stream.readline()
    if pdb_fasta == "":
        logger.error("The chain %s has not been found in %s", chain, pdbid)
        raise RuntimeError
    return pdb_fasta


def _get_chain(line):
    if not line.find(":"):
        logger.error("No chain information in %s", line)
        raise RuntimeError
    return line[6]


def print_database_uniprots(database):
    from libs.data_model import UniprotEntry
    i = 0
    print("*" * 80)
    for entry in database.iter(UniprotEntry):
        print(entry)
        i += 1
    print(i, "<- Database size")
    print("*" * 80)


CACHE_DIRECTORY = "./cache"


class UrlCache(object):
    """
    Cache HTTP GET request.
    TODO Add remove need to specify name and enable use of arbitrary URL.
    """
    
    def __init__(self, identifier, url, cache_name="url-cache", mode="r"):
        self.identifier = str(identifier)
        self.url = url
        self.cache_directory = os.path.join(CACHE_DIRECTORY, cache_name)
        os.makedirs(self.cache_directory, exist_ok=True)
        self.mode = mode
        self.stream = None
    
    def __enter__(self):
        if not os.path.exists(self._path()):
            self.download()
        self.stream = open(self._path(), self.mode)
        return self.stream
    
    def __exit__(self, *args):
        self.stream.close()
    
    def _path(self):
        return os.path.join(self.cache_directory, self.identifier)
    
    def as_path(self):
        if not os.path.exists(self._path()):
            self.download()
        return os.path.join(self.cache_directory, self.identifier)
    
    def download(self):
        with urllib.request.urlopen(self.url) as response:
            with open(self._path(), "wb") as output_stream:
                output_stream.write(response.read())

