#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path

base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir_bigbwt = ""
data_dir_pfp    = ""

data_sizes = [64,128,256]
w_values   = [10,15,20,25,30]
p_values   = [50,100,150,200]
f_values   = [0.001,0.005,0.01,0.05,0.1,0.5]
n_threads  = 32

#------------------------------------------------------------
# execute command: return True if everything OK, False otherwise
def execute_command(command, seconds):
    try:
        print("Executing: {}".format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid)
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        return False
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        print("Command exceeded timeout:")
        print("\t"+ command)
        return False
    return True


def get_pscan(out_dir):
    repository = "https://github.com/alshai/Big-BWT.git"
    get_command = "git clone {} {}".format(repository, out_dir + "/Big-BWT")
    execute_command(get_command, 100)
    build_command = "make -C {}".format(out_dir + "/Big-BWT")
    execute_command(build_command, 1000)


def main():
    get_pscan(base_dir)
    pscan_exe = base_dir + "/Big-BWT/newscan.x"
    pfp_exe = "pfp++"

    # Run bigbwt
    base_command = "/usr/bin/time --verbose {pscan} -w {c_w} -p {c_p} -t {c_threads} -f {file}"
    for size in data_sizes:
        for w in w_values:
            for p in p_values:
                command = base_command.format(pscan = pscan_exe, c_w = w, c_p = p, c_threads = n_threads,
                                              file = data_dir_bigbwt + "/chr19.{}.fa".format(size))
                execute_command(command, 1000000)

    # Run pfp
    base_command = "/usr/bin/time --verbose {pfp} -v {vcf} -r {ref} -o {out} -m {c_size} -w {c_w} -p {c_p} -f {c_f} -t {c_threads} -s --use-acceleration --print-statistics"
    for size in data_sizes:
        for w in w_values:
            for p in p_values:
                for f in f_values:
                    command = base_command.format(
                        pfp = pfp_exe,
                        vcf = data_dir_pfp + "/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
                        ref = data_dir_pfp + "/Homo_sapiens.GRCh37.dna.chromosome.19.fa.gz",
                        out = data_dir_pfp + "/out_w{}_p{}_f{}_s{}".format(w,p,f,size),
                        c_size = size, c_w = w, c_p = p, c_f = f, c_threads = n_threads)
                    execute_command(command, 1000000)


if __name__ == '__main__':
    main()