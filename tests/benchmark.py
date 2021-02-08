#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, wget

base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir_bigbwt = ""
data_dir_pfp    = ""

DEBUG = True
if (DEBUG):
    data_sizes = [32]
    w_values   = [10]
    p_values   = [50]
    f_values   = [0.1]
    n_threads  = 8
else:
    data_sizes = [64,128,256,512,1000, 2000]
    w_values   = [10,15,20,25,30]
    p_values   = [50,100,150,200]
    f_values   = [0.001,0.005,0.01,0.05,0.1,0.5]
    n_threads  = 32

#------------------------------------------------------------
# execute command: return command's stdout if everything OK, None otherwise
def execute_command(command, seconds):
    try:
        print("Executing: {}".format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid, stdout=subprocess.PIPE)
        (output, err) = process.communicate()
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        return None
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        print("Command exceeded timeout:")
        print("\t"+ command)
        return None
    return output


def get_pscan(out_dir):
    repository = "https://github.com/alshai/Big-BWT.git"
    get_command = "git clone {} {}".format(repository, out_dir + "/Big-BWT")
    execute_command(get_command, 100)
    build_command = "make -C {}".format(out_dir + "/Big-BWT")
    execute_command(build_command, 1000)

def get_vcf_file(out_dir, chromosome_id):
    if chromosome_id not in [str(c) for c in range(1,23)]:
        print('Chromosome id currently not supported')
        return
    chromosome_name = "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502." \
                       "genotypes.vcf.gz".format(chromosome_id)
    chromosome_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/" \
                     + chromosome_name
    print('Downloading {}'.format(chromosome_url))
    wget.download(chromosome_url, out_dir + '/' + chromosome_name)

def extract_fasta(out_dir, vcf_file, n_sequences=1):
    vcf_to_fasta_exe = "vcf_to_fa"


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
                out = execute_command(command, 1000000)

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
                    out = execute_command(command, 1000000)


if __name__ == '__main__':
    main()