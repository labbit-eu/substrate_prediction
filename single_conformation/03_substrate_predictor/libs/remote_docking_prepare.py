import os
import logging
import shutil
import tarfile
import stat

from libs import paths
from libs import docking
from libs.utils import remove_existing
from libs.utils import read_json
from libs.configuration import ENV, REMOTE, TASK
from external import remote_processing as remote

logger = logging.getLogger(__name__)


def prepare_docking_batches():
    run_docking_path = os.path.join(paths.temp(), "run_docking")
    if os.path.exists(run_docking_path):
        shutil.rmtree(run_docking_path)
    remove_existing(os.path.join(paths.temp(), "VS_storage.json"))

    files_in_batches, max_docking_cost = _assign_files_in_subdirs_to_batches(paths.candidates_pdbqt())

    logger.info("\n***************************************************************************************************"
                "*********\nAn estimated execution time for the most costly candidate batch is %ds per CPU core."
                "\nMake sure that the queue and walltime specified in your template file allows for completion "
                "including\ngenerous error margin as the execution time is extrapolated from the "
                "substrate docking on the local host.\n***************************************************************"
                "*********************************************", round(max_docking_cost*int(TASK["docking"]["cpu"])))
    file_template = create_docking_batch_script("vina_template.txt")
    processed_files = 0

    for batch_id in range(0, len(files_in_batches)):
        batch_name = "batch{:=04d}".format(batch_id+1)
        batch_dir = os.path.join(run_docking_path, batch_name)

        os.makedirs(batch_dir)
        os.makedirs(os.path.join(batch_dir, "log"))
        os.makedirs(os.path.join(batch_dir, "out"))
        os.makedirs(os.path.join(batch_dir, "in"))
        shutil.copy(paths.target_protein(), os.path.join(run_docking_path, batch_dir, "in", "protein.pdbqt"))
        shutil.copy(paths.candidates_vina_template(),
                    os.path.join(run_docking_path, batch_dir, "in", "vina_template.txt"))

        file_path = os.path.join(run_docking_path, batch_dir, "run_docking.sh")
        with open(file_path, "w") as output_stream:
            output_stream.write(file_template.replace("$$REMOTE.BATCH_NAME$$", batch_name))

        os.chmod(file_path, stat.S_IRWXU)

        if REMOTE["selfsustained"]:
            shutil.copy(ENV["vina"], os.path.join(run_docking_path, batch_dir, "in", "vina"))
        else:
            if REMOTE["vina_dir"] is None:
                message = "When parameter \"selfsustained\" in [REMOTE] section of the configuration file is set to " \
                          "False, parameter \"vina_dir\" in [REMOTE] has to be defined and has to contain location " \
                          "of vina binary on the remote host."
                raise RuntimeError(message)
            else:
                remote_connection = remote.RemoteHandler()
                vina_file_info = remote_connection.file_info(os.path.join(REMOTE["vina_dir"], "vina"))
                if vina_file_info is None:
                    if remote_connection.make_folder_path(REMOTE["vina_dir"]):
                        remote_connection.send_file(ENV["vina"], os.path.join(REMOTE["vina_dir"], "vina"), mode=0o500)
                    else:
                        raise RuntimeError("Cannot create {} folder since an object of same name exists on remote host."
                                           " Try changing \"vina_dir\" parameter in [REMOTE] section of the"
                                           " configuration file.".format(REMOTE["vina_dir"]))

        batch_tar_file = tarfile.open(os.path.join(run_docking_path, batch_name, "pdbqt.tar.gz"), "w:gz")
        for file_path in files_in_batches[batch_id]:
            subdir_name, file_name = paths.get_subdir_filename_from_filepath(file_path)
            os.makedirs(os.path.join(batch_dir, "out", subdir_name), exist_ok=True)
            os.makedirs(os.path.join(batch_dir, "log", subdir_name), exist_ok=True)
            batch_tar_file.add(file_path, os.path.join("pdbqt", subdir_name, file_name))

        batch_tar_file.close()
        processed_files += len(files_in_batches[batch_id])
    logger.info("In total, %d files were processed.", processed_files)


def _get_batch_command(config_filename):
    return """for subdir in `ls pdbqt`; do
  for f in `ls pdbqt/$subdir/*.pdbqt`; do
    b=`basename $f .pdbqt`
    echo $b
    {} --config in/{} --ligand $f --out out/$subdir/${{b}}.pdbqt --log log/$subdir/${{b}}.log
  done
done  
""".format(os.path.join(REMOTE["vina_dir"], "vina"), config_filename)


def _get_batch_template():
    if REMOTE["queue"]["template_file"] is not None:
        logger.info("Using the user defined template file %s for remote scripts.", REMOTE["queue"]["template_file"])
        with open(REMOTE["queue"]["template_file"], "r") as input_stream:
            return input_stream.read()
    else:
        if REMOTE["queue"]["manager"] == "slurm":
            logger.info("Using the default SLURM template file for remote scripts.")
            return """#!$$REMOTE.SHELL_PATH$$
#SBATCH --nodes=1
#SBATCH --cpus-per-task=$$REMOTE.CPU$$
#SBATCH --partition=$$REMOTE.QUEUE_NAME$$
#SBATCH --mem-per-cpu=$$REMOTE.MEM$$
#SBATCH --job-name=$$REMOTE.BATCH_NAME$$

# set-up working directory
export TMPDIR=$$REMOTE.TMPDIR$$
mkdir -p ${TMPDIR}
cd $TMPDIR

# copy data to working directory
cp -r ${SLURM_SUBMIT_DIR}/* ${TMPDIR}
tar -xf pdbqt.tar.gz
rm pdbqt.tar.gz

# execute screening
$$autogenerated_code$$

# copy data to submit directory
tar -czf out.tar.gz out
tar -czf log.tar.gz log
cp -r $TMPDIR/*.tar.gz $SLURM_SUBMIT_DIR

# cleaning - this is deactivated on purpose, review and activate or do yourself
#rm -rf $TMPDIR
"""
        elif REMOTE["queue"]["manager"] == "torque":
            return """#!$$REMOTE.SHELL_PATH$$
#PBS -l nodes=1:ppn=$$REMOTE.CPU$$
#PBS -q $$REMOTE.QUEUE_NAME$$
#PBS -l mem=$$REMOTE.MEM$$
#PBS -N $$REMOTE.BATCH_NAME$$

# set-up working directory
export TMPDIR=$$REMOTE.TMPDIR$$
mkdir -p ${TMPDIR}
cd $TMPDIR

# copy data to working directory
cp -r ${PBS_O_WORKDIR}/* ${TMPDIR}
tar -xf pdbqt.tar.gz
rm pdbqt.tar.gz

# execute screening
$$autogenerated_code$$

# copy data to submit directory
tar -czf out.tar.gz out
tar -czf log.tar.gz log
cp -r $TMPDIR/*.tar.gz $PBS_O_WORKDIR

# cleaning - this is deactivated on purpose, review and activate or do yourself
#rm -rf $TMPDIR
"""


def create_docking_batch_script(config_filename):
    template = _get_batch_template().replace("$$REMOTE.SHELL_PATH$$", REMOTE["queue"]["shell_path"])
    template = template.replace("$$REMOTE.CPU$$", str(TASK["docking"]["cpu"]))
    template = template.replace("$$REMOTE.QUEUE_NAME$$", REMOTE["queue"]["name"])
    if REMOTE["queue"]["manager"] == "slurm":
        template = template.replace("$$REMOTE.MEM$$", str(REMOTE["queue"]["mem"]))
    elif REMOTE["queue"]["manager"] == "torque":
        template = template.replace("$$REMOTE.MEM$$", str(REMOTE["queue"]["mem"] * TASK["docking"]["cpu"]))
    template = template.replace("$$REMOTE.TMPDIR$$", REMOTE["queue"]["tmpdir"])
    template = template.replace("$$autogenerated_code$$", _get_batch_command(config_filename))

    return template


def _convert_batch_size_number(num_files, batch_size_or_number):
    if num_files % batch_size_or_number == 0:
        return num_files // batch_size_or_number
    else:
        return num_files // batch_size_or_number + 1


def _get_batch_number(number_pdbqt_files):
    if REMOTE["batch_size"] is not None:
        batch_size = REMOTE["batch_size"]
        if batch_size < 1:
            raise RuntimeError("Failed to prepare batches of size {} for remote execution. Parameter \"batch_size\" in "
                               "[REMOTE] section of the configuration file must be integer >= 1\n".format(batch_size))
        batch_number = _convert_batch_size_number(number_pdbqt_files, batch_size)
    else:
        batch_number = int(REMOTE["batch_number"])
        if batch_number < 1:
            raise RuntimeError("Failed to prepare {} batches for virtual screening of candidates on remote host. "
                               "Parameter \"batch_number\" in [REMOTE] section of the configuration file must be "
                               "integer >= 1\n".format(batch_number))
        batch_size = _convert_batch_size_number(number_pdbqt_files, batch_number)

    if batch_size < 50:
        logger.warning(
            "Currently creating %d batches of size %d.\n Too small batch sizes are frequently rather inefficient, "
            "consider increasing the  \"batch_size\" or decreasing \"batch_number\" in section REMOTE.",
            batch_number, batch_size)

    logger.info("Splitting %d candidates into %d batches for remote docking", number_pdbqt_files, batch_number)
    return batch_number


def _is_resolved(subdir_name, file_name, resolved_candidates):
    for batch_name in resolved_candidates.keys():
        if subdir_name not in resolved_candidates[batch_name].keys():
            continue
        if file_name in resolved_candidates[batch_name][subdir_name].keys():
            return True

    return False


def _assign_files_in_subdirs_to_batches(in_dir):
    files = {}
    num_files = 0
    resolved_candidates = {}
    cost_function = read_json(paths.candidates_info())["cost_function"]

    if os.path.exists(os.path.join(paths.candidates(), "resolved_failed_candidates.json")):
        resolved_candidates = read_json(os.path.join(paths.candidates(), "resolved_failed_candidates.json"))

    for subdir_name in os.listdir(in_dir):
        subdir = os.path.join(in_dir, subdir_name)
        for file_name in os.listdir(subdir):
            if os.path.exists(os.path.join(paths.candidates_docking(), subdir_name, file_name)) or \
                    _is_resolved(subdir_name, file_name, resolved_candidates):
                continue

            pdbqt_file = os.path.join(subdir, file_name)
            cost = cost_function[docking.get_torsdof_from_pdbqt(pdbqt_file)]
            if cost not in files.keys():
                files[cost] = []
            files[cost].append(pdbqt_file)
            num_files += 1

    if num_files == 0:
        return [], 0
    else:
        return _split_files_by_docking_cost(files, _get_batch_number(num_files))


def _get_file_with_maxcost(files, costs):
    for cost in costs:
        if len(files[cost]):
            return files[cost].pop(), cost

    return None


def _split_files_by_docking_cost(files, batch_number):
    costs = sorted(files.keys(), reverse=True)
    cost_of_batches = {}
    files_in_batches = {}

    for i in range(1, batch_number+1):
        cost_of_batches[i] = 0
        files_in_batches[i] = []

    for cost in costs:
        num_files = len(files[cost])
        for i in range(num_files // batch_number):
            for j in cost_of_batches.keys():
                files_in_batches[j].append(files[cost].pop())
                cost_of_batches[j] += cost

    file_with_maxcost = True

    while file_with_maxcost:
        for j in range(1, batch_number+1):
            file_with_maxcost = _get_file_with_maxcost(files, costs)
            if file_with_maxcost is not None:
                files_in_batches[j].append(file_with_maxcost[0])
                cost_of_batches[j] += file_with_maxcost[1]
            else:
                break
        if file_with_maxcost is not None:
            for j in range(batch_number, 0, -1):
                file_with_maxcost = _get_file_with_maxcost(files, costs)
                if file_with_maxcost is not None:
                    files_in_batches[j].append(file_with_maxcost[0])
                    cost_of_batches[j] += file_with_maxcost[1]
                else:
                    break

    files_to_return = []
    cost_of_batches_to_return = []

    for i in range(1, batch_number+1):
        files_to_return.append(files_in_batches[i])
        cost_of_batches_to_return.append(cost_of_batches[i])

    return files_to_return, max(cost_of_batches_to_return)
