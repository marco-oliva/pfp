#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, wget

base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir_bigbwt = "/blue/boucher/marco.oliva/tmp"
#data_dir_pfp    = "/blue/boucher/marco.oliva/data/1kgp"
data_dir_pfp    = "/blue/boucher/marco.oliva/tmp/generated_vcf"

# VCF list
vcf_list = ''
#for i in range(1,23):
#    vcf_list += '"' + data_dir_pfp + '/vcf/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",'.format(i)
#vcf_list += '"' + data_dir_pfp + '/vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz",'
#vcf_list += '"' + data_dir_pfp + '/vcf/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz",'
#vcf_list += '"' + data_dir_pfp + '/vcf/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"'
#vcf_list += '"' + data_dir_pfp + '/vcf/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"'.format(19)
vcf_list += '"' + data_dir_pfp + '/generated.vcf.gz"'

# Reference list
ref_list = ''
#for i in range(1,23):
#    ref_list += '"' + data_dir_pfp + '/reference/{}.fa.gz",'.format(i)
#ref_list += '"' + data_dir_pfp + '/reference/X.fa.gz",'
#ref_list += '"' + data_dir_pfp + '/reference/Y.fa.gz",'
#ref_list += '"' + data_dir_pfp + '/reference/MT.fa.gz"'
ref_list += '"' + data_dir_pfp + '/22.fa.gz"'

# Generate config file
with open(base_dir + '/config.ini', 'w') as config_file:
    config_file.write('# Generated config file\n')
    config_file.write('vcf = [{}]\n'.format(vcf_list))
    config_file.write('ref = [{}]\n'.format(ref_list))


DEBUG = False
RUN_BIGBWT = False
RUN_WITH_MODULO = False
RUN_WITHOUT_MODULO = True
if (DEBUG):
    data_sizes = [1]
    w_values   = [10]
    p_values   = [50]
    f_values   = [1.0]
    n_threads  = 8
else:
    data_sizes = [500]
    w_values   = [10, 20, 30, 40]
    p_values   = [100, 500, 1000, 2000]
    f_values   = [1.0, 0.01]
    F_values   = [0.5]
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

def extract_fasta(out_dir, ref_list, vcf_list, n_sequences=1):
    vcf_to_fasta_exe = "../build/vcf_to_fa"
    command = "/usr/bin/time --verbose {} --configure {} -m {} -o {}".format(vcf_to_fasta_exe, base_dir + '/config.ini',
                                                                             n_sequences, out_dir + '/extracted.fa')
    execute_command(command, 100000000)

def main():
    pfp_exe = base_dir + "/../build/pfp++"

    # Run bigbwt
    if RUN_BIGBWT:
        get_pscan(data_dir_bigbwt)
        pscan_exe = data_dir_bigbwt + "/Big-BWT/newscan.x"
        base_command = "/usr/bin/time --verbose {pscan} -w {c_w} -p {c_p} -t {c_threads} -f {file}"
        for size in data_sizes:
            extract_fasta(data_dir_bigbwt, ref_list, vcf_list, size)
            for w in w_values:
                for p in p_values:
                    command = base_command.format(pscan = pscan_exe, c_w = w, c_p = p, c_threads = n_threads,
                                                  file = data_dir_bigbwt + "/extracted.fa")
                    out = execute_command(command, 1000000)

    # Run pfp, with modulo
    if (RUN_WITH_MODULO):
        base_command = "/usr/bin/time --verbose {pfp} --configure {conf} -o {out} -m {c_size} -w {c_w} -p {c_p} -f {c_f} -F {c_F} -t {c_threads} --seeds --use-acceleration --print-statistics"
        for size in data_sizes:
            for w in w_values:
                for p in p_values:
                    for f in f_values:
                        for F in F_values:
                            command = base_command.format(
                                pfp = pfp_exe,
                                conf = base_dir + '/config.ini',
                                out = base_dir + "/out_w{}_p{}_f{}_F{}_s{}".format(w,p,f,F,size),
                                c_size = size, c_w = w, c_p = p, c_f = f, c_F = F, c_threads = n_threads)
                            out = execute_command(command, 1000000)
                            with open('out_w{}_p{}_f{}_F{}_s{}.log'.format(w,p,f,F,size), 'wb') as log_file:
                                log_file.write(out)

    # Run pfp, without modulo
    if (RUN_WITHOUT_MODULO):
        base_command = "/usr/bin/time --verbose {pfp} --configure {conf} -o {out} -m {c_size} -w {c_w} --not-use-modulo -f {c_f} -F {c_F} -t {c_threads} --seeds --use-acceleration --print-statistics"
        for size in data_sizes:
            for w in w_values:
                for f in f_values:
                    for F in F_values:
                        command = base_command.format(
                            pfp = pfp_exe,
                            conf = base_dir + '/config.ini',
                            out = base_dir + "/out_w{}_f{}_F{}_s{}".format(w,f,F,size),
                            c_size = size, c_w = w, c_f = f, c_F = F, c_threads = n_threads)
                        out = execute_command(command, 1000000)
                        with open('out_w{}_f{}_F{}_s{}.log'.format(w,f,F,size), 'wb') as log_file:
                            log_file.write(out)

if __name__ == '__main__':
    main()
