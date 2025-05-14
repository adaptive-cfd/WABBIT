#!/usr/bin/env python3
import os, sys, argparse, subprocess, time, shutil, glob, logging, select

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# these are all tests that we can run
# new groups can be added with "---NAME---" and new tests by providing name, wavelet and dimension
# new tests need to be implemented in the class WabbitTest by doing:
#   - defining the test folder location in __init__
#   - defining what should be executed in run
#   - defining where the log-file should be in init_logging

group_names = ["post", "wavelets", "ghost_nodes", "invertibility", "adaptive", "convection", "acm", "insects", "cvs"]
tests = [
        # f"---{group_names[0]}---",
        # "TESTING/wabbit_post/pod/pod_test.sh",
        f"---{group_names[1]}---",  # group identifier
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF20", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF22", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF24", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF26", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF28", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF40", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF42", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF44", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF46", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF60", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF62", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF64", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF66", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF20", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF22", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF40", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF42", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF44", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF60", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF62", "dim":3},

        f"---{group_names[2]}---",  # group identifier
        {"test_name":"ghost_nodes", "wavelet":"CDF20", "dim":2},
        {"test_name":"ghost_nodes", "wavelet":"CDF40", "dim":2},
        {"test_name":"ghost_nodes", "wavelet":"CDF60", "dim":2},
        {"test_name":"ghost_nodes", "wavelet":"CDF20", "dim":3},
        {"test_name":"ghost_nodes", "wavelet":"CDF40", "dim":3},
        {"test_name":"ghost_nodes", "wavelet":"CDF60", "dim":3},

        f"---{group_names[3]}---",  # group identifier
        {"test_name":"invertibility", "wavelet":"CDF20", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF22", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF24", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF26", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF28", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF40", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF42", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF44", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF46", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF60", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF62", "dim":2},
        {"test_name":"invertibility", "wavelet":"CDF64", "dim":2},

        f"---{group_names[4]}---",  # group identifier
        {"test_name":"adaptive", "wavelet":"CDF20", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF22", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF40", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF42", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF44", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF60", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF62", "dim":2},

        f"---{group_names[5]}---",  # group identifier
        {"test_name":"blob_equi", "wavelet":"CDF40", "dim":2},
        {"test_name":"blob_equi", "wavelet":"CDF20", "dim":3},
        {"test_name":"blob_equi", "wavelet":"CDF40", "dim":3},
        {"test_name":"blob_adaptive", "wavelet":"CDF20", "dim":2},
        {"test_name":"blob_adaptive", "wavelet":"CDF22", "dim":2},
        {"test_name":"blob_adaptive", "wavelet":"CDF40", "dim":2},
        {"test_name":"blob_adaptive", "wavelet":"CDF42", "dim":2},
        {"test_name":"blob_adaptive", "wavelet":"CDF44", "dim":2},
        {"test_name":"blob_adaptive", "wavelet":"CDF60", "dim":2},
        {"test_name":"blob_adaptive", "wavelet":"CDF62", "dim":2},
        {"test_name":"blob_adaptive", "wavelet":"CDF22", "dim":3},
        {"test_name":"blob_adaptive", "wavelet":"CDF40", "dim":3},
        {"test_name":"blob_adaptive", "wavelet":"CDF44", "dim":3},

        f"---{group_names[6]}---",  # group identifier
        {"test_name":"acm", "wavelet":"CDF40", "dim":2},
        {"test_name":"acm", "wavelet":"CDF44", "dim":2},
        {"test_name":"acm_norm", "wavelet":"CDF44", "dim":2},
        {"test_name":"acm_significant", "wavelet":"CDF44", "dim":2},

        f"---{group_names[7]}---",  # group identifier
        {"test_name":"dry_fractal_tree", "wavelet":"CDF22", "dim":3},
        {"test_name":"dry_bumblebee", "wavelet":"CDF22", "dim":3},
        {"test_name":"dry_emundus_4wings", "wavelet":"CDF22", "dim":3},
        {"test_name":"dry_muscaComplete", "wavelet":"CDF22", "dim":3},
        {"test_name":"dry_dipteraFourier", "wavelet":"CDF22", "dim":3},
        {"test_name":"dry_dipteraHermite", "wavelet":"CDF22", "dim":3},
        {"test_name":"dry_dipteraBodyRotation", "wavelet":"CDF22", "dim":3},
        {"test_name":"dry_paratuposaComplete", "wavelet":"CDF22", "dim":3},

#        f"---{group_names[8]}---",  # group identifier
#        {"test_name":"denoise_butterfly", "wavelet":"CDF42", "dim":2},
#        {"test_name":"denoise_grey", "wavelet":"CDF42", "dim":2},
    ]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


fail_color, pass_color, end_color, underline = '\033[31;1m', '\033[92;1m', '\033[0m', '\033[4m'
import importlib  
try: wabbit_tools = importlib.import_module("wabbit_tools")  # this needs file with wabbit_tools loaded to env PYTHONPATH
except:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Since 15 Aug 2023, the unit testing framework has evolved. It now stores full HDF5 files in the TESTING directory, which makes it easier to visualize the reference data and current results, should they be different.")
    print("We now calculate the L2 error of the field, if the grid is identical. This new framework requires the https://github.com/adaptive-cfd/python-tools repository for comparing two WABBIT HDF5 files.")
    print(f"\n You do not seem to have {fail_color}wabbit_tools{end_color} available! Either you do not have the repository, or its directory is not in your $PYTHONPATH")
    print(f"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(f"{fail_color}Cannot run unit tests !!{end_color}")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    sys.exit(881)

def run_command(command, logger):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8', errors='ignore')

    # Poll both stdout and stderr
    while True:
        reads = [process.stdout.fileno(), process.stderr.fileno()]
        ret = select.select(reads, [], [])

        for fd in ret[0]:
            if fd == process.stdout.fileno():
                read = process.stdout.readline()
                if read: logger.info(read.strip("\n"))
            if fd == process.stderr.fileno():
                read = process.stderr.readline()
                if read: logger.error(read.strip('\00').strip("\n"))
        
        # very little wait so that line buffers can be filled appropriately, elsewise
        time.sleep(0.01)
                    
        # Check if the process has finished
        if process.poll() is not None: break

    # Ensure all remaining output is processed
    for read in process.stdout:
        if read: logger.info(read.strip("\n"))

    for read in process.stderr:
        if read: logger.error(read.strip('\00').strip("\n"))
    
    return process

# this object defines a wabbit test
class WabbitTest:
    ini = test_name = wavelet = dim = mpi_command = memory = test_dir = cwd = run_dir = logger = None
    valid = True  # check if something in initialization went wrong

    # init the class itself
    def __init__(self, ini=None, test_name=None, wavelet=None, dim=None, mpi_command="nice mpirun -n 4", memory="8.0GB", run_dir=None):
        if run_dir == None: self.run_dir = os.getcwd()  # this should be the run directory
        else: self.run_dir = os.path.abspath(run_dir)
        self.ini = ini
        self.test_name = test_name
        self.wavelet = wavelet
        self.dim = dim
        self.mpi_command = mpi_command
        self.memory = memory

        # define here where the test folder will be located relative to the run directory!
        if self.test_name in ["equi_refineCoarsen_FWT_IWT", "ghost_nodes", "invertibility"]:
            self.test_dir = os.path.join(self.run_dir, "TESTING", "wavelets")
        elif self.test_name == "adaptive":
            self.test_dir = os.path.join(self.run_dir, "TESTING", "wavelets", f"{self.test_name}_{self.wavelet}")
        elif self.test_name in ["blob_equi", "blob_adaptive"]:
            self.test_dir = os.path.join(self.run_dir, "TESTING", "conv", f"{self.test_name}_{self.dim}D_{self.wavelet}")
        elif self.test_name in ["acm", "acm_norm", "acm_significant"]:
            self.test_dir = os.path.join(self.run_dir, "TESTING", "acm", f"{self.test_name}_{self.wavelet}")
        elif self.test_name in ["dry_fractal_tree", "dry_muscaComplete", "dry_dipteraFourier", "dry_dipteraHermite", "dry_dipteraBodyRotation", "dry_bumblebee", "dry_emundus_4wings", "dry_paratuposaComplete"]:
            self.test_dir = os.path.join(self.run_dir, "TESTING", "insects", f"{self.test_name}_{self.wavelet}")
        elif self.test_name in ["denoise_butterfly", "denoise_grey"]:
            self.test_dir = os.path.join(self.run_dir, "TESTING", "cvs", f"{self.test_name}_{self.wavelet}")
        else:
            print("Unknown test name")
            self.valid = False
        if not os.path.isdir(self.test_dir):
            print(f"Skipping a test. Test directory does not exist: {self.test_dir}")
            self.valid = False

        if self.valid:
            os.chdir(self.test_dir)


    # actually run the test! Here new tests can be added, when working with ini-files, use the block with blob and acm and point to the correct ini file
    def run(self, write_diff=False):
        if self.test_name == "equi_refineCoarsen_FWT_IWT":
            # change to directory
            command1 = f"{self.mpi_command} {self.run_dir}/wabbit-post --refine-coarsen-test --wavelet={self.wavelet} --memory={self.memory} --dim={self.dim}"
            command2 = f"{self.mpi_command} {self.run_dir}/wabbit-post --wavelet-decomposition-unit-test --wavelet={self.wavelet} --memory={self.memory} --dim={self.dim}"
            result1 = run_command(command1, self.logger)
            if result1.returncode != 0:
                return result1
            result2 = run_command(command2, self.logger)
            return result2
        
        elif self.test_name == "ghost_nodes":
            command = f"{self.mpi_command} {self.run_dir}/wabbit-post --ghost-nodes-test --wavelet={self.wavelet} --memory={self.memory} --dim={self.dim}"
            result = run_command(command, self.logger)
            return result
        
        elif self.test_name == "invertibility":
            command = f"{self.mpi_command} {self.run_dir}/wabbit-post --wavelet-decomposition-invertibility-test --wavelet={self.wavelet} --memory={self.memory} --dim={self.dim}"
            result = run_command(command, self.logger)
            return result
        
        elif self.test_name == "adaptive":
            in_file = os.path.join("..", "..", "vor_000020000000.h5")  # relative to tmp_dir
            
            # change to directory to tmp
            tmp_dir = f"{self.test_dir}/tmp"
            if not os.path.isdir(tmp_dir): os.mkdir(tmp_dir)
            os.chdir(tmp_dir)

            # run commands
            file1, file2 = "vor_00100.h5", "vor_00200.h5"
            command1 = f"{self.mpi_command} {self.run_dir}/wabbit-post --refine-everywhere \"{in_file}\" \"{file1}\" --wavelet={self.wavelet} --time=1.0"
            command2 = f"{self.mpi_command} {self.run_dir}/wabbit-post --coarsen-everywhere \"{file1}\" \"{file2}\" --wavelet={self.wavelet} --time=2.0"
            result1 = run_command(command1, self.logger)
            if result1.returncode != 0:
                return result1
            result2 = run_command(command2, self.logger)

            # the refined file is quite data-heavy but always the same betweend different Y-wavelets if CDFXY. Let's just keep them and test them for unlifted wavelets!
            if (self.wavelet not in ["CDF20", "CDF40", "CDF60"]):
                os.remove("vor_00100.h5")

            # compare all files present in test_dir
            try:
                all_similar = self.compare_files(tmp_dir, write_diff=write_diff)
                result2.returncode = not all_similar
            # catch any error - for example HDF5 error - and then say the test failed
            # this is important to still print the log (and continue)
            except:
                result2.returncode = 1

            # change back to test_dir
            os.chdir(self.test_dir)
            return result2
        elif "denoise" in self.test_name:
            # change to directory to tmp
            tmp_dir = f"{self.test_dir}/tmp"
            if not os.path.isdir(tmp_dir): os.mkdir(tmp_dir)
            os.chdir(tmp_dir)

            test_file = self.test_name.split("_")[1]

            # define commands, the first creates the file to be denoised and the second actually does the denoising
            denoise_file = f"../{test_file}.png"
            denoise_h5 = f"./{test_file}.h5"
            noise_add = "-n 10"if test_file == "grey" else ""  # grey file gets noise added on top
            command1 = f"image2wabbit.py {denoise_file} -o {denoise_h5} --level 5 --bs 16 {noise_add}"  # image is already noisy so no extra noise is added
            command2 = f"{self.mpi_command} {self.run_dir}/wabbit-post --denoise --files=\"{denoise_h5}\" --wavelet={self.wavelet} --memory={self.memory}"
                        
            # first, convert image to a valid wabbit file
            result1 = run_command(command1, self.logger)
            if result1.returncode != 0:
                return result1
            # now, denoise the file
            result2 = run_command(command2, self.logger)

            # # remove files which are not used for comparisons
            # os.remove(denoise_h5)

            # compare all files present in test_dir
            try:
                all_similar = self.compare_files(tmp_dir, write_diff=write_diff)
                result1.returncode = not all_similar
            # catch any error - for example HDF5 error - and then say the test failed
            # this is important to still print the log (and continue)
            except Exception as e:
                self.logger.error("ERROR: Was not able to compare files")
                self.logger.error(e)
                result1.returncode = 1

            # change back to test_dir
            os.chdir(self.test_dir)
            return result1
        # this part is meant for any tests which simply call an ini file, just provide the ini-file in the beginning and the rest is handled automatically
        elif self.test_name in ["dry_fractal_tree", "dry_muscaComplete", "dry_dipteraFourier", "dry_dipteraHermite", "dry_dipteraBodyRotation", "dry_bumblebee", "dry_emundus_4wings", "dry_paratuposaComplete"]:
            ini_file = os.path.join("..", "PARAMS_dry_run.ini")  # relative to tmp_dir

            # change to directory to tmp
            tmp_dir = f"{self.test_dir}/tmp"
            if not os.path.isdir(tmp_dir): os.mkdir(tmp_dir)
            os.chdir(tmp_dir)

            save_us = ""
            if self.test_name in ["dry_emundus_4wings", "dry_muscaComplete", "dry_dipteraFourier", "dry_dipteraHermite", "dry_dipteraBodyRotation", "dry_paratuposaComplete"]:
                save_us = "--save-us"

            # run simmulation
            command1 = f"{self.mpi_command} {self.run_dir}/wabbit-post --dry-run {ini_file} --memory={self.memory} --pruned {save_us}"
            result1 = run_command(command1, self.logger)

            # compare all files present in test_dir
            try:
                all_similar = self.compare_files(tmp_dir, write_diff=write_diff)
                result1.returncode = not all_similar
            # catch any error - for example HDF5 error - and then say the test failed
            # this is important to still print the log (and continue)
            except Exception as e:
                self.logger.error("ERROR: Was not able to compare files")
                self.logger.error(e)
                result1.returncode = 1

            # change back to test_dir
            os.chdir(self.test_dir)
            return result1
        elif self.test_name in ["blob_equi", "blob_adaptive", "acm", "acm_norm", "acm_significant"]:
            # lets say where the ini-file is
            if self.test_name in ["blob_equi", "blob_adaptive"]:
                ini_file = os.path.join("..", "blob-convection.ini")  # relative to tmp_dir
            elif self.test_name in ["acm", "acm_norm", "acm_significant"]:
                ini_file = os.path.join("..", "acm_cyl.ini")  # relative to tmp_dir

            # change to directory to tmp
            tmp_dir = f"{self.test_dir}/tmp"
            if not os.path.isdir(tmp_dir): os.mkdir(tmp_dir)
            os.chdir(tmp_dir)

            # run simmulation
            command1 = f"{self.mpi_command} {self.run_dir}/wabbit {ini_file} --memory={self.memory}"
            result1 = run_command(command1, self.logger)

            # compare all files present in test_dir
            try:
                all_similar = self.compare_files(tmp_dir, write_diff=write_diff)
                result1.returncode = not all_similar
            # catch any error - for example HDF5 error - and then say the test failed
            # this is important to still print the log (and continue)
            except Exception as e:
                self.logger.error("ERROR: Was not able to compare files")
                self.logger.error(e)
                result1.returncode = 1

            # change back to test_dir
            os.chdir(self.test_dir)
            return result1
        else:
            self.logger.error("Not implemented yet")
            return False
    

    # I want to log at the same time to the console and possibly files as well, so I solve this with the logging module which handles the streams
    # for a new test the log-file location needs to be specified
    def init_logging(self, verbose=False, suite_log_handler=None, stdout_handler=None):
        self.log_file = None
        if self.test_name in ["equi_refineCoarsen_FWT_IWT", "ghost_nodes", "invertibility"]:
            self.log_file = os.path.join(self.test_dir, f"{self.test_name}_{self.dim}D_{self.wavelet}.log")
        elif self.test_name == "adaptive":
            self.log_file = os.path.join(self.test_dir, "run.log")
        elif self.test_name in ["blob_equi", "blob_adaptive"]:
            self.log_file = os.path.join(self.test_dir, "blob-convection.log")
        elif self.test_name == ["acm", "acm_norm", "acm_significant"]:
            self.log_file = os.path.join(self.test_dir, "acm_cyl.log")
        elif self.test_name in ["dry_fractal_tree", "dry_muscaComplete", "dry_dipteraFourier", "dry_dipteraHermite", "dry_bumblebee", "dry_emundus_4wings", "dry_dipteraBodyRotation", "dry_paratuposaComplete"]:
            self.log_file = os.path.join(self.test_dir, "dry_run.log")
        elif self.test_name in ["denoise_butterfly", "denoise_grey"]:
            self.log_file = os.path.join(self.test_dir, "denoise.log")
        else:
            self.log_file = os.path.join(self.test_dir, "test.log")
        
        open(self.log_file, 'w').close()  # Clear the log file at the beginning
        # Set up logging for runs, which writes to different log file and with verbose to test log file and stdout as well
        self.logger = logging.getLogger("logger_run")
        self.logger.setLevel(logging.INFO)
        test_log_handler = logging.FileHandler(self.log_file, mode='a')
        test_log_handler.setFormatter(logging.Formatter('%(message)s'))
        self.logger.handlers = []
        self.logger.addHandler(test_log_handler)
        if verbose:
            self.logger.addHandler(suite_log_handler)
            self.logger.addHandler(stdout_handler)


    # remove residual files and possibly overwrite reference results
    def clean_up(self, replace=False, keep_tmp=False, logger=logger):
        if self.test_name in ["equi_refineCoarsen_FWT_IWT", "ghost_nodes", "invertibility"]:
            # remove files - only .dat files are created
            if not keep_tmp:
                for file in glob.glob("*.dat"):
                    os.remove(file)
        else:
            if replace:
                # remove all old files
                for file in glob.glob(os.path.join(self.test_dir, "*.h5")):
                    os.remove(file)
                    logger.info(f"   Del  ref file: {file}")
                # copy over new files, we are currently in tmp folder
                for file in glob.glob(os.path.join(self.test_dir, "tmp", "*.h5")):
                    shutil.copy(file, os.path.join(self.test_dir, os.path.split(file)[1]))
                    logger.info(f"   Copy new file: {file}")
                    logger.info(f"         to file: {os.path.join(self.test_dir, os.path.split(file)[1])}")
            if not keep_tmp:
                shutil.rmtree(os.path.join(self.test_dir, "tmp"))
        
        # change back to directory where we were
        os.chdir(self.run_dir)


    # takes every reference file in the test folder and tries to compare it to available test results in tmp folder
    def compare_files(self, tmp_dir, verbose=True, write_diff=False):
        tmp_files = glob.glob(os.path.join(tmp_dir, "*.h5"))
        tmp_split = [os.path.split(i_file)[1] for i_file in tmp_files]
        all_similar = True
        if verbose:
            self.logger.info("\n" * 10 + "=" * 80 + "\nRun done, analyzing results now\n" + "=" * 80 + "\n" * 10)
        happy, sad = 0, 0  # count the number of tests
        for file in sorted(glob.glob(os.path.join(self.test_dir, "*.h5"))):
                file_split = os.path.split(file)[1]
                if file_split in tmp_split:
                    t_compare_s = time.time()
                    if verbose:
                        self.logger.info("*" * 80 + "\nComparing wabbit HDF5 files")
                        self.logger.info(f"Reference file (1) = {os.path.join(self.test_dir, file_split)}")
                        self.logger.info(f"Test result    (2) = {os.path.join(tmp_dir, file_split)}")

                    state_ref = wabbit_tools.WabbitHDF5file()
                    state_ref.read(file, verbose=False)
                    state_new = wabbit_tools.WabbitHDF5file()
                    state_new.read(os.path.join(tmp_dir, file_split), verbose=False)

                    bool_similar = state_ref.isClose(state_new, verbose=verbose, logger=self.logger)
                    if not bool_similar:
                        all_similar = False
                        sad += 1
                    else: happy += 1

                    # this is to be able to visualize the difference
                    if not bool_similar and write_diff:
                        state_diff = state_ref - state_new
                        state_diff.write(os.path.join(self.test_dir, "diff-" + file_split))
                        # state_new_int = (state_ref * 0) + state_new  # sneaky way to interpolate w_obj2 grid onto w_obj_new
                        # state_new_int.write(os.path.join(self.test_dir, "new-" + file_split))
                    
                    if verbose:
                        t_compare = time.time() - t_compare_s 
                        self.logger.info(f"Finished comparison. Time= {t_compare:7.3f} s\n")
                # file is not present in tmp directory
                else:
                    # it could start with "new" or "diff" - this is likely a debug file used for visualization so lets skip it
                    if file_split.startswith("new") or file_split.startswith("diff"): continue

                    self.logger.info("*" * 80 + "\nComparing wabbit HDF5 files")
                    self.logger.info(f"Reference file (1) = {os.path.join(self.test_dir, file_split)}")
                    self.logger.info(f"ERROR: file not found for test result")
                    all_similar = False
                    sad += 1

        # write summary
        self.logger.info(f"\nFinished {happy+sad} tests")
        self.logger.info(f"\tHappy tests: {happy:2d} :)")
        self.logger.info(f"\t  Sad tests: {sad:2d} :(")
        return all_similar



def main():
    parser = argparse.ArgumentParser(description='Run WABBIT unit tests.')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print output of test directly to screen')
    mpi_group = parser.add_mutually_exclusive_group()
    mpi_group.add_argument('--nprocs', type=int, default=4, help='Number of processors, default is 4')
    mpi_group.add_argument('--mpi_command', type=str, default=None, help=r'MPI command, default is "nice mpirun -n {nprocs}"')
    parser.add_argument('--memory', type=str, default='8.0GB', help='Memory flag, default is \"8.0GB\"')
    parser.add_argument('--test', type=str, default="all", help="Specific test to run, provide \"NAME-DIM-WAVELET\". You can also provide \"all\" or a specific group")
    # parser.add_argument('--print_test', action='store_true', help="Print output of test directly to screen")
    rep_group = parser.add_mutually_exclusive_group()
    rep_group.add_argument('--replace', action='store_true', help="Replace reference values with results")
    rep_group.add_argument('--replace-fail', action='store_true', help="Replace reference values with results only for failing tests")
    parser.add_argument('--keep-tmp', action="store_true", help="Do not delete the temporary directories to investigate files")
    parser.add_argument('--write-diff', action="store_true", help="Write the difference between reference and deviating new results to a file. Could be possibly interpolated")
    parser.add_argument('--wabbit-dir', type=str, default=None, help="Location to where WABBIT is located to run the tests if not in present directory")
    args = parser.parse_args()

    nprocs = args.nprocs
    mpi_command = args.mpi_command or f"nice mpirun -n {nprocs}"
    memory = args.memory

    # check run directory
    if args.wabbit_dir == None: args.wabbit_dir = os.getcwd()  # this should be the run directory
    else: args.wabbit_dir = os.path.abspath(args.wabbit_dir)
    # check if we can access wabbit-post
    executable = shutil.which(os.path.join(args.wabbit_dir, "wabbit"))
    if executable is None:
        logging.error(f"ERROR: Did not find wabbit executable in wabbit_dir: {args.wabbit_dir}")
        sys.exit(1)
    
    # Set up logging for general output
    # this at the same time writes to the log file and to stdout
    log_total = os.path.join(args.wabbit_dir, "TESTING", "test.log")
    open(log_total, 'w').close()  # Clear the log file at the beginning
    logger_suite = logging.getLogger("logger_suite")
    logger_suite.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    suite_log_handler = logging.FileHandler(log_total, mode='a')
    suite_log_handler.setFormatter(formatter)
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(formatter)
    logger_suite.addHandler(suite_log_handler)
    logger_suite.addHandler(stdout_handler)

    # check if user actually wants to replace test results
    if args.replace or args.replace_fail:
        response = input("Are you REALLY sure you want to reset all or some reference results? (yes/no): ").strip().lower()
        if response in ['yes', 'y']:
            logger_suite.info("Going to replace all or some reference results!")
        else:
            logger_suite.info("Ok alright, I understand. I also often do not know what I am doing so come back once you know that it's the right time.")
            sys.exit()

    # run all tests
    group_run = "all"  # init group variable
    if args.test == "all":
        tests_run = tests
        logger_suite.info(f"\n\t{underline}WABBIT: run all existing unit tests{end_color}\n")
    # run a group of tests only
    elif args.test in group_names:
        group_run = args.test
        tests_run = tests  # we will check later what tests to run from this group
    else:
        # check if this test exists
        test_parts = args.test.split("-")
        if len(test_parts) != 3:
            test_exists = False
        # replace D in dimension entry if it was provided
        try:
            test_parts[1] = test_parts[1].replace("D", "")
            test_dict = {"test_name":test_parts[0], "wavelet": test_parts[2], "dim": int(test_parts[1])}
            test_exists = test_dict in tests
        except: test_exists = False

        if test_exists:
            logger_suite.info(f"\n\t \033[4m WABBIT: run existing unit test {args.test} \033[0m\n")
            tests_run = [test_dict]
        else:
            logger_suite.error(f"ERROR: {args.test} is not a valid unit test or group of tests")
            sys.exit(1)

    happy_sum = 0
    sad_sum = 0
    summary = []

    # give user some information
    logger_suite.info(f"employed command for parallel exec: {mpi_command}")
    logger_suite.info(f"memory flag for wabbit is: {memory}\n")

    start_time = time.time()

    group_now = ""
    for ts in tests_run:
        if isinstance(ts, str):
            if ts.startswith("---"):
                group_now = ts.replace("-", "")
                if group_run == "all" or group_now == group_run: print(ts)
        else:
            # cycle if this is not part of the group we want to check
            if group_run != "all" and group_now != group_run: continue

            # build test object
            test_obj = WabbitTest(test_name=ts["test_name"], dim=ts["dim"], wavelet=ts["wavelet"], mpi_command=mpi_command, memory=memory, run_dir=args.wabbit_dir)
            if test_obj.valid == False: continue  # in case initialization did not work out

            # create logger for this specific test
            test_obj.init_logging(verbose=args.verbose, suite_log_handler=suite_log_handler, stdout_handler=stdout_handler)
            # print some infos to console about this test, disable line-break so that it looks nice
            for i_handler in logger_suite.handlers: i_handler.terminator = ""
            logger_suite.info(f"Test= {ts['test_name']:<50}wavelet={ts['wavelet']} dim={ts['dim']}D")
            for i_handler in logger_suite.handlers: i_handler.terminator = "\n"

            ts_start_time = time.time()
            result = test_obj.run(write_diff=args.write_diff)
            ts_end_time = time.time() - ts_start_time

            for i_handler in logger_suite.handlers: i_handler.terminator = ""
            if isinstance(result, subprocess.Popen):
                # write output to console, make it a bit fancy
                if result.returncode == 0:
                    print(f"{pass_color}", end="")  # this only works for console
                    logger_suite.info(f"\tPass ")
                    happy_sum += 1
                    summary.append(0)
                else:
                    print(f"{fail_color}", end="")  # this only works for console
                    logger_suite.info(f"\tFail ")
                    sad_sum += 1
                    summary.append(1)

                for i_handler in logger_suite.handlers: i_handler.terminator = "\n"
                print(f"{end_color}", end="")  # this only works for console                    
                logger_suite.info(f"\tTime= {ts_end_time:7.3f} s")

                if result.returncode != 0:
                    # test failed, provide direct link to log file
                    print(f"{end_color}", end="")  # this only works for console  
                    logger_suite.info(f"logfile: \t"+test_obj.log_file+"\n")
            else:
                print(f"{fail_color}", end="")  # this only works for console
                logger_suite.info(f"\tFail ")
                print(f"{end_color}", end="")  # this only works for console   

                for i_handler in logger_suite.handlers: i_handler.terminator = "\n"
                logger_suite.info(f"\tTest was not executed")

            # remove temporary dir and replace reference results if wanted
            if args.replace_fail:
                test_obj.clean_up(replace=(result.returncode != 0), keep_tmp=args.keep_tmp, logger=logger_suite)
            else:
                test_obj.clean_up(replace=args.replace, keep_tmp=args.keep_tmp, logger=logger_suite)
            

    total_time = time.time() - start_time
    logger_suite.info(f"\nFinished all tests. Time= {total_time:7.3f} s\n")

    # give a little summary
    logger_suite.info(f"All in all we have:")
    for i_handler in logger_suite.handlers: i_handler.terminator = ""
    logger_suite.info(f"\t")
    for i_res in summary:
        if i_res==0:
            print(f"{pass_color}", end="")  # this only works for console
            logger_suite.info(f"O")
        else:
            print(f"{fail_color}", end="")  # this only works for console
            logger_suite.info(f"X")
        print(f"{end_color}", end="")  # this only works for console                    
    for i_handler in logger_suite.handlers: i_handler.terminator = "\n"

    print(f"{pass_color}", end="")  # this only works for console
    logger_suite.info(f"\n\n\t   sum happy tests:\t{happy_sum}")
    print(f"{fail_color}", end="")  # this only works for console
    logger_suite.error(f"\t   sum sad tests:\t{sad_sum}")
    print(f"{end_color}", end="")  # this only works for console                    


if __name__ == "__main__":
    main()
