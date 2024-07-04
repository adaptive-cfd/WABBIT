#!/usr/bin/env python3
import numpy as np
import os, sys, argparse, subprocess, select, logging, shutil
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": "true",
    "font.size": 11,
    "axes.linewidth": 1,
    "lines.linewidth": 1,
    'figure.dpi': 141
})

wavelet_list = ["CDF20", "CDF22", "CDF40", "CDF42", "CDF44", "CDF60", "CDF62"]
# wavelet_list = ["CDF44"]
log_file = "output.log"
# Set up logging
logging.basicConfig(level=logging.INFO, format='%(message)s',
                    handlers=[
                        logging.FileHandler(log_file, mode='a'),
                        logging.StreamHandler(sys.stdout)
                    ])

def run_tests(args):
    for wavelet in wavelet_list:
        process = subprocess.Popen(args.mpi_command.split(" ") + [os.path.join(args.wabbit_dir, "wabbit-post"), "--compression-unit-test",
            f"--wavelet={wavelet}", f"--Jmax={args.Jmax}", f"--dim={args.dim}", "--eps-norm=Linfty", f"--memory={args.memory}", "--eps-normalized=1", f"--Bs={args.Bs}"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        error = False

        # Poll both stdout and stderr
        while True:
            reads = [process.stdout.fileno(), process.stderr.fileno()]
            ret = select.select(reads, [], [])

            for fd in ret[0]:
                if fd == process.stdout.fileno():
                    read = process.stdout.readline()
                    if read and args.verbose: logging.info(read.strip())
                    elif any(substring in read for substring in ["eps=", "Elapsed", "Nb="]):
                        logging.info(read.strip())
                if fd == process.stderr.fileno():
                    read = process.stderr.readline()
                    if read:
                        logging.error(read.rstrip('\00').strip())
                        error = True
                        
            # Check if the process has finished
            if process.poll() is not None:
                break

        # Ensure all remaining output is processed
        for read in process.stdout:
            if read and args.verbose: logging.info(read.strip())
            elif any(substring in read for substring in ["eps=", "Elapsed"]):
                logging.info(read.strip())

        for read in process.stderr:
            if read:
                logging.error(read.rstrip('\00').strip())
                error = True

        # move 
        if error == False:
            shutil.move("error.csv", f"error_{wavelet}.csv")


def plot(args):
    fig1 = plt.figure(1, figsize=[4, 4])
    fig2 = plt.figure(2, figsize=[4, 4])
    fig3 = plt.figure(3, figsize=[4, 4])

    for wavelet in wavelet_list:  
        d = np.loadtxt(f'error_{wavelet}.csv', delimiter=';', skiprows=1)
        
        fig1.gca().loglog( d[:,0], d[:,1], '.-' , label =f"$L_2$ {wavelet}", mfc='w', markersize=5)
        fig2.gca().loglog( d[:,0], d[:,2], '.-' , label =f"$L_\infty$ {wavelet}", mfc='w', markersize=5)
        
        ####  plot compression rate
        fig3.gca().loglog(d[:,0], d[:,3] / 2**(args.dim*args.Jmax), '.-', label = wavelet, mfc='w', markersize=5)

    fig1.gca().loglog(d[:,0], d[:,0], 'k--')
    fig2.gca().loglog(d[:,0], d[:,0], 'k--')
    fig1.gca().set_ylabel("Compression error")
    fig2.gca().set_ylabel("Compression error")
    fig3.gca().set_ylabel(f"Compression rate to {2**(args.dim*args.Jmax)} blocks")
    for i_fig in [fig1, fig2, fig3]:
        plt.figure(i_fig)
        plt.legend(handlelength=1, columnspacing=1, labelspacing=0.1, handletextpad=0.3)
        plt.xlabel("Threshold epsilon")
        plt.suptitle(f"Bs=${args.Bs}$ Jmax=${args.Jmax}$ dim=${args.dim}$ eps-norm=$L_\infty$ eps-normalized=T", fontsize=10)
        plt.tight_layout(rect=[0, 0.01, 1, 0.98], pad=0.15)  # rect respects suptitle

    plt.figure(fig1); plt.savefig("compression-error-L2.png")
    plt.figure(fig2); plt.savefig("compression-error-Linfty.png")
    plt.figure(fig3); plt.savefig("compression-ratio.png")
    try:
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages('compression-test.pdf') as pdf:
            plt.figure(fig1); pdf.savefig()
            plt.figure(fig2); pdf.savefig()
            plt.figure(fig3); pdf.savefig()
    except Exception as e:
        print("An error occurred when trying to creatre multi-page pdf document:")
        print(e)


def main():
    parser = argparse.ArgumentParser(description='Run compression test.')
    mpi_group = parser.add_mutually_exclusive_group()
    mpi_group.add_argument('--nprocs', type=int, default=4, help='Number of processors, default is 4')
    mpi_group.add_argument('--mpi_command', type=str, default=None, help=r'MPI command, default is "nice mpirun -n {nprocs}"')
    parser.add_argument('--memory', type=str, default='8.0GB', help='Memory flag, default is \"8.0GB\"')
    parser.add_argument("--Jmax", type=int, default=6, help="Maximum depth, gives (2**dim)**Jmax blocks, defaults to 6.")
    parser.add_argument("-d", "--dim", type=int, default=2, help="Dimension, 2 or 3, defaults to 2.")
    parser.add_argument("--Bs", type=int, default=24, help="Blocksize, gives Bs*Bs points per block, defaults to 24.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Output full log file")
    parser.add_argument("-p", "--plot", action="store_true", help="Plot results")
    parser.add_argument("-s", "--skip-run", action="store_true", help="Skip execution of tests")
    parser.add_argument('--wabbit-dir', type=str, default=None, help="Location to where WABBIT is located to run the tests if not in present directory")
    args = parser.parse_args()

    args.mpi_command = args.mpi_command or f"nice mpirun -n {args.nprocs}"
    memory = args.memory
    if args.wabbit_dir == None: args.wabbit_dir = os.getcwd()  # this should be the run directory
    else: args.wabbit_dir = os.path.abspath(args.wabbit_dir)

    # Clear the log file at the beginning
    open(log_file, 'w').close()

    logging.info(f"employed command for parallel exec: {args.mpi_command}")
    logging.info(f"memory flag for wabbit is: {args.memory}\n")

    # check if we can access wabbit-post
    executable = shutil.which(os.path.join(args.wabbit_dir, "wabbit-post"))
    if executable is None:
        logging.error(f"ERROR: Did not find wabbit executable in wabbit_dir: {args.wabbit_dir}")
        sys.exit(1)

    if not args.skip_run:
        run_tests(args)

    if args.plot:
        plot(args)
    
    os.remove("header.txt")  # this is not really needed for this version

if __name__ == "__main__":
    main()