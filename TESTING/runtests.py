#!/usr/bin/env python3
import os, sys, argparse, subprocess, time, shutil, glob, io


fail_color, pass_color, end_color = '\033[31;1m', '\033[92;1m', '\033[0m'
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

# this object defines a wabbit test
class WabbitTest:
    ini = test_name = wavelet = dim = mpi_command = memory = test_dir = cwd = run_dir = None
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

        if self.test_name in ["equi_refineCoarsen_FWT_IWT", "ghost_nodes"]:
            self.test_dir = os.path.join(self.run_dir, "TESTING", "wavelets")
        elif self.test_name == "adaptive":
            self.test_dir = os.path.join(self.run_dir, "TESTING", "wavelets", f"{self.test_name}_{self.wavelet}")
        elif self.test_name in ["blob_equi", "blob_adaptive"]:
            self.test_dir = os.path.join(self.run_dir, "TESTING", "conv", f"{self.test_name}_{self.dim}D_{self.wavelet}")
        elif self.test_name == "acm":
            self.test_dir = os.path.join(self.run_dir, "TESTING", "acm", f"acm_cyl_adaptive_{self.wavelet}")
        else:
            print("Unknown test name")
            self.valid = False
        if not os.path.isdir(self.test_dir):
            print(f"Skipping a test. Test directory does not exist: {self.test_dir}")
            self.valid = False

        if self.valid:
            os.chdir(self.test_dir)

    def run(self):
        cwd = os.getcwd()  # this should be the run directory
        if self.test_name == "equi_refineCoarsen_FWT_IWT":
            # change to directory
            command1 = f"{self.mpi_command} {self.run_dir}/wabbit-post --refine-coarsen-test --wavelet={self.wavelet} --memory={self.memory} --dim={self.dim}"
            command2 = f"{self.mpi_command} {self.run_dir}/wabbit-post --wavelet-decomposition-unit-test --wavelet={self.wavelet} --memory={self.memory} --dim={self.dim}"
            result1 = subprocess.run(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            if result1.returncode != 0:
                return result1
            result2 = subprocess.run(command2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            # combine the logs
            result2.stdout = result1.stdout + "\n" + result2.stdout
            return result2
        elif self.test_name == "ghost_nodes":
            command = f"{self.mpi_command} {self.run_dir}/wabbit-post --ghost-nodes-test --wavelet={self.wavelet} --memory={self.memory} --dim={self.dim}"
            result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            return result
        elif self.test_name == "adaptive":
            in_file = os.path.join("..", "..", "vor_000020000000.h5")  # relative to tmp_dir
            
            # change to directory to tmp
            tmp_dir = f"{self.test_dir}/tmp"
            if not os.path.isdir(tmp_dir): os.mkdir(tmp_dir)
            os.chdir(tmp_dir)

            # run commands
            file1, file2 = "vor_00100.h5", "vor_00200.h5"
            command1 = f"{self.mpi_command} {self.run_dir}/wabbit-post --sparse-to-dense \"{in_file}\" \"{file1}\" --wavelet={self.wavelet} --operator=refine-everywhere"
            command2 = f"{self.mpi_command} {self.run_dir}/wabbit-post --sparse-to-dense \"{file1}\" \"{file2}\" --wavelet={self.wavelet} --operator=coarsen-everywhere"
            result1 = subprocess.run(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            if result1.returncode != 0:
                return result1
            result2 = subprocess.run(command2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            # combine the logs
            result2.stdout = result1.stdout + "\n" + result2.stdout

            # compare all files present in test_dir
            try:
                all_similar, output = self.compare_files(tmp_dir)
                result2.stdout += output
                result2.returncode = not all_similar
            # catch any error - for example HDF5 error - and then say the test failed
            # this is important to still print the log (and continue)
            except:
                result2.returncode = 1

            # change back to test_dir
            os.chdir(self.test_dir)
            return result2
        elif self.test_name in ["blob_equi", "blob_adaptive", "acm"]:
            if self.test_name in ["blob_equi", "blob_adaptive"]:
                ini_file = os.path.join("..", "blob-convection.ini")  # relative to tmp_dir
            elif self.test_name in ["acm"]:
                ini_file = os.path.join("..", "acm_cyl.ini")  # relative to tmp_dir

            # change to directory to tmp
            tmp_dir = f"{self.test_dir}/tmp"
            if not os.path.isdir(tmp_dir): os.mkdir(tmp_dir)
            os.chdir(tmp_dir)

            # run simmulation
            command1 = f"{self.mpi_command} {self.run_dir}/wabbit {ini_file} --memory={self.memory}"
            result1 = subprocess.run(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')

            # compare all files present in test_dir
            try:
                all_similar, output = self.compare_files(tmp_dir)
                result1.stdout += output
                result1.returncode = not all_similar
            # catch any error - for example HDF5 error - and then say the test failed
            # this is important to still print the log (and continue)
            except:
                result1.returncode = 1

            # change back to test_dir
            os.chdir(self.test_dir)
            return result1
        else:
            print("Not implemented yet")
            return False
    
    def get_log_file(self):
        if self.test_name in ["equi_refineCoarsen_FWT_IWT", "ghost_nodes"]:
            return os.path.join(self.test_dir, f"{self.test_name}_{self.dim}D_{self.wavelet}.log")
        elif self.test_name == "adaptive":
            return os.path.join(self.test_dir, "run.log")
        elif self.test_name in ["blob_equi", "blob_adaptive"]:
            return os.path.join(self.test_dir, "blob-convection.log")
        elif self.test_name == "acm":
            return os.path.join(self.test_dir, "acm_cyl.log")


    def clean_up(self, replace=False):
        if self.test_name in ["equi_refineCoarsen_FWT_IWT", "ghost_nodes"]:
            # remove files - only .dat files are created
            for file in glob.glob("*.dat"):
                os.remove(file)
        else:
            if replace:
                # remove all old files
                for file in glob.glob(os.path.join(self.test_dir, "*.h5")):
                    os.remove(file)
                    print(f"   Del  ref file: {file}")
                # copy over new files, we are currently in tmp folder
                for file in glob.glob(os.path.join(self.test_dir, "tmp", "*.h5")):
                    shutil.copy(file, os.path.join(self.test_dir, os.path.split(file)[1]))
                    print(f"   Copy new file: {file}")
                    print(f"         to file: {os.path.join(self.test_dir, os.path.split(file)[1])}")

            shutil.rmtree(os.path.join(self.test_dir, "tmp"))
        
        # change back to directory where we were
        os.chdir(self.run_dir)

    def compare_files(self, tmp_dir, verbose=True, write_diff=False):
        tmp_files = glob.glob(os.path.join(tmp_dir, "*.h5"))
        tmp_split = [os.path.split(i_file)[1] for i_file in tmp_files]
        all_similar = True
        # redirect stdout to variable as we want to put infos into logfile
        output = io.StringIO()
        sys.stdout = output
        if verbose:
            print("\n" * 10 + "=" * 80 + "\nRun done, analyzing results now\n" + "=" * 80 + "\n" * 10)
        happy, sad = 0, 0  # count the number of tests
        for file in sorted(glob.glob(os.path.join(self.test_dir, "*.h5"))):
                file_split = os.path.split(file)[1]
                if file_split in tmp_split:
                    t_compare_s = time.time()
                    if verbose:
                        print("*" * 80 + "\nComparing wabbit HDF5 files")
                        print(f"Reference file (1) = {os.path.join(self.test_dir, file_split)}")
                        print(f"Test result    (2) = {os.path.join(tmp_dir, file_split)}")

                    state_ref = wabbit_tools.WabbitHDF5file()
                    state_ref.read(file, verbose=False)
                    state_new = wabbit_tools.WabbitHDF5file()
                    state_new.read(os.path.join(tmp_dir, file_split), verbose=False)

                    bool_similar = state_ref.isClose(state_new, verbose=verbose)
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
                        print(f"Finished comparison. Time= {t_compare:7.3f} s\n")
                # file is not present in tmp directory
                else:
                    # it could start with "new" or "diff" - this is likely a debug file used for visualization so lets skip it
                    if file_split.startswith("new") or file_split.startswith("diff"): continue

                    print("*" * 80 + "\nComparing wabbit HDF5 files")
                    print(f"Reference file (1) = {os.path.join(self.test_dir, file_split)}")
                    print(f"ERROR: file not found for test result")
                    all_similar = False
                    sad += 1

        # write summary
        print(f"\nFinished {happy+sad} tests")
        print(f"\tHappy tests: {happy:2d} :)")
        print(f"\t  Sad tests: {sad:2d} :(")

        # reset stdout
        sys.stdout = sys.__stdout__
        return all_similar, output.getvalue()
                    

# these are all tests that we can run
tests = [
        # "---post---",
        # "TESTING/wabbit_post/pod/pod_test.sh",
        "---wavelets---",  # group identifier
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF20", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF20", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF22", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF40", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF42", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF44", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF60", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF62", "dim":2},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF20", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF22", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF40", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF42", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF44", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF60", "dim":3},
        {"test_name":"equi_refineCoarsen_FWT_IWT", "wavelet":"CDF62", "dim":3},
        {"test_name":"ghost_nodes", "wavelet":"CDF20", "dim":2},
        {"test_name":"ghost_nodes", "wavelet":"CDF40", "dim":2},
        {"test_name":"ghost_nodes", "wavelet":"CDF60", "dim":2},
        {"test_name":"ghost_nodes", "wavelet":"CDF20", "dim":3},
        {"test_name":"ghost_nodes", "wavelet":"CDF40", "dim":3},
        {"test_name":"ghost_nodes", "wavelet":"CDF60", "dim":3},
        {"test_name":"adaptive", "wavelet":"CDF20", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF22", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF40", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF42", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF44", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF60", "dim":2},
        {"test_name":"adaptive", "wavelet":"CDF62", "dim":2},

        "---convection---",  # group identifier
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

        "---acm---",  # group identifier
        {"test_name":"acm", "wavelet":"CDF44", "dim":2},
    ]



def main():
    parser = argparse.ArgumentParser(description='Run WABBIT unit tests.')
    mpi_group = parser.add_mutually_exclusive_group()
    mpi_group.add_argument('--nprocs', type=int, default=4, help='Number of processors, default is 4')
    mpi_group.add_argument('--mpi_command', type=str, default=None, help=r'MPI command, default is "nice mpirun -n {nprocs}"')
    parser.add_argument('--memory', type=str, default='8.0GB', help='Memory flag, default is \"8.0GB\"')
    parser.add_argument('--test', type=str, default="all", help="Specific test to run, provide \"NAME-DIM-WAVELET\". You can also provide \"all\" or a specific group")
    parser.add_argument('--print_test', action='store_true', help="Print output of test directly to screen")
    rep_group = parser.add_mutually_exclusive_group()
    rep_group.add_argument('--replace', action='store_true', help="Replace reference values with results")
    rep_group.add_argument('--replace-fail', action='store_true', help="Replace reference values with results only for failing tests")
    parser.add_argument('--wabbit-dir', type=str, default=None, help="Location to where WABBIT is located to run the tests if not in present directory")
    args = parser.parse_args()

    nprocs = args.nprocs
    mpi_command = args.mpi_command or f"nice mpirun -n {nprocs}"
    memory = args.memory

    if args.replace:
        response = input("Are you REALLY sure you want to reset all test results? (yes/no): ").strip().lower()
        if response in ['yes', 'y']:
            print("Going to replace all reference results!")
        else:
            print("Ok alright, I understand. I also often do not know what I am doing so come back once you know that it's the right time.")
            sys.exit()

    # run all tests
    group_run = "all"  # init group variable
    if args.test == "all":
        tests_run = tests
        print("\n\t \033[4m WABBIT: run all existing unit tests \033[0m\n")
    # run a group of tests only
    elif args.test in ["post", "wavelets", "convection", "acm"]:
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
            print(f"\n\t \033[4m WABBIT: run existing unit test {args.test} \033[0m\n")
            tests_run = [test_dict]
        else:
            print(f"ERROR: {args.test} is not a valid unit test")
            sys.exit(1)

    happy_sum = 0
    sad_sum = 0
    summary = []

    print(f"employed command for parallel exec: {mpi_command}")
    print(f"memory flag for wabbit is: {memory}\n")
    # ToDo!
    # print("to modify the command, pass --memory=[MEMORY] and/or --mpi_command=[MPI COMMAND] in shell\n")

    if nprocs != 4:
        print(f"{fail_color}WARNING{end_color}")
        print("your tests might fail because the keyvalues for load balancing may differ if you don't use nprocs=4 for testing")

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

            logfile = test_obj.get_log_file()
            if os.path.exists(logfile): os.remove(logfile)

            print(f"Test= {ts['test_name']:<50}wavelet={ts['wavelet']} dim={ts['dim']}D", end="", flush=True)

            ts_start_time = time.time()
            result = test_obj.run()
            ts_end_time = time.time() - ts_start_time

            if isinstance(result, subprocess.CompletedProcess):
                with open(logfile, 'w') as f:
                    f.write(result.stdout)
                    f.write(result.stderr)
                if args.print_test:
                    print(result.stdout)
                    print(result.stderr)

                if result.returncode == 0:
                    print(f"{pass_color}\tPass {end_color}\tTime= {ts_end_time:7.3f} s")
                    happy_sum += 1
                    summary.append(0)
                else:
                    print(f"{fail_color}\tFail {end_color}\tTime= {ts_end_time:7.3f} s")
                    sad_sum += 1
                    summary.append(1)
            else:
                print(f"{fail_color}\tFail {end_color}\tTest was not executed")

            if args.replace_fail:
                test_obj.clean_up(replace=(result.returncode != 0))
            else:
                test_obj.clean_up(replace=args.replace)
            

    total_time = time.time() - start_time
    print(f"\nFinished all tests. Time= {total_time:7.3f} s\n")
    # print("All in all we have:\n")

    # for i, ts in enumerate(tests):
    #     if not ts.startswith("---"):
    #         status = f"{pass_color}ok{end_color}" if summary[i] == 0 else f"{fail_color}X{end_color}"
    #         print(f"{ts:80s} {status}")

    # print(f"\n\t {pass_color}sum happy tests: {end_color}\t{happy_sum}")
    # print(f"\n\t {fail_color}sum sad tests: {end_color}\t{sad_sum}")

if __name__ == "__main__":
    main()
