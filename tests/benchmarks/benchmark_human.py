#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, wget, errno, datetime, logging
from Bio import SeqIO

# Set this dir
data_base_dir = '/home/marco/tmp/tests'
pre_download_data_dir = '/home/marco/Data/1kgp'
version_letter = 'a'

# Assuming installed
# - bcftools
# - htslib
# - gzip
# - bgzip

project_base_dir = os.path.dirname(os.path.abspath(__file__))
date_string             = datetime.datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
data_dir_bigbwt         = "{}/human_tests/{}/bigbwt".format(data_base_dir, date_string)
data_dir_pfp            = "{}/human_tests/{}/pfp".format(data_base_dir, date_string)
common_data_dir         = "{}/human_tests/data".format(data_base_dir)
common_tools_dir        = "{}/human_tests/tools".format(data_base_dir)
tmp_fasta_dir           = "{}/human_tests/data/tmp".format(data_base_dir)

# Generate config file
# with open(project_base_dir + '/config.ini', 'w') as config_file:
#     config_file.write('# Generated config file\n')
#     config_file.write('vcf = [{}]\n'.format(vcf_list))
#     config_file.write('ref = [{}]\n'.format(ref_list))


data_sizes = [10, 100, 200]
w_values   = [10]
p_values   = [100]
n_threads  = 32

#------------------------------------------------------------
# execute command: return command's stdout if everything OK, None otherwise
def execute_command(command, time_it=False, seconds=1000000):
    try:
        if time_it:
            command = '/usr/bin/time --verbose {}'.format(command)
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


def get_pscan(work_dir):
    repository = "https://github.com/alshai/Big-BWT.git"
    get_command = "git clone {} {}".format(repository, work_dir + "/Big-BWT")
    execute_command(get_command)
    build_command = "make -C {}".format(work_dir + "/Big-BWT")
    execute_command(build_command)
    return work_dir + '/Big-BWT/pscan.x'

def get_au_pair(work_dir):
    return "binary_path"

def get_pfp(work_dir):
    repository = "https://github.com/marco-oliva/pfp.git"
    execute_command("git clone {} {}".format(repository, work_dir + "/pfp"))
    mkdir_p(work_dir + '/pfp/build')
    os.chdir(work_dir + '/pfp/build')
    execute_command("cmake .. ")
    execute_command("make -j")
    return work_dir + '/pfp/build/pfp++'

def get_vcf_files(out_dir):
    vcf_files_list = list()
    for chromosome_id in [str(c) for c in range(1,23)]:
        chromosome_name = "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5{}.20130502." \
                          "genotypes.vcf.gz".format(chromosome_id, version_letter)
        print(pre_download_data_dir + '/vcf/' + chromosome_name)
        if (os.path.exists(pre_download_data_dir + '/vcf/' + chromosome_name)):
            print('Copying {}'.format(chromosome_name))
            execute_command('cp {} {}'.format(pre_download_data_dir + '/vcf/' + chromosome_name,
                                              out_dir + '/' + chromosome_name)
                            )
        else:
            chromosome_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/" \
                             + chromosome_name
            print('Downloading {}'.format(chromosome_url))
            wget.download(chromosome_url, out_dir + '/' + chromosome_name)

        print('From gzip to bgzip')
        execute_command('gunzip {}'.format(out_dir + '/' + chromosome_name), time_it=True)
        execute_command('bgzip {}'.format(out_dir + '/' + chromosome_name[-3]), time_it=True)
        execute_command('bcftools index {}'.format(out_dir + '/' + chromosome_name[-3] + '.bgz'), time_it=True)
        vcf_files_list.append(out_dir + '/' + chromosome_name[-3] + '.bgz')
    return vcf_files_list

def get_reference_files(out_dir):
    reference_url = 'ftp://ftp.1000genomes.ebi.ac.uk//vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
    if (os.path.exists(pre_download_data_dir + '/reference/hs37d5.fa.gz')):
        execute_command('cp {} {}'.format(pre_download_data_dir + '/reference/hs37d5.fa.gz', out_dir + 'hs37d5.fa.gz'))
    else:
        print("Downloading reference files")
        wget.download(reference_url, out_dir + '/' + 'hs37d5.fa.gz')
    print('Extracting chromosomes')
    execute_command("gunzip {}".format(out_dir + '/' + 'hs37d5.fa.gz'))

    out_fasta_list = list()
    for record in SeqIO.parse(out_dir + '/' + 'hs37d5.fa', 'fasta'):
        out_fasta_list.append(out_dir + '/' + record.id)
        with gzip.open(out_dir + '/' + record.id, 'wb') as output_fasta:
            SeqIO.write(record, output_fasta, "fasta")
    return out_fasta_list

def extract_fasta(out_file_path, ref_file, vcf_file, sample_id):
    bcf_command = "/usr/bin/time --verbose bcftools consensus -H 1 -f {} -s {} {} > {}".format(
        red_file, sample_id, vcf_file, out_file_path)
    execute_command(bcf_command, time_it=True)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python ≥ 2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise # nop

def main_debug():
    vcf_dir = common_data_dir + '/vcf'
    mkdir_p(vcf_dir)
    print(get_vcf_files(vcf_dir))


# def main():
#     pfp_exe = project_base_dir + "/../build/pfp++"
#
#     # Run bigbwt
#     if RUN_BIGBWT:
#         get_pscan(data_dir_bigbwt)
#         pscan_exe = data_dir_bigbwt + "/Big-BWT/newscan.x"
#         base_command = "/usr/bin/time --verbose {pscan} -w {c_w} -p {c_p} -t {c_threads} -f {file}"
#         for size in data_sizes:
#             extract_fasta(data_dir_bigbwt, ref_list, vcf_list, size)
#             for w in w_values:
#                 for p in p_values:
#                     command = base_command.format(pscan = pscan_exe, c_w = w, c_p = p, c_threads = n_threads,
#                                                   file = data_dir_bigbwt + "/extracted.fa")
#                     out = execute_command(command, 1000000)
#
#     # Run pfp, with modulo
#     if (RUN_WITH_MODULO):
#         base_command = "/usr/bin/time --verbose {pfp} --configure {conf} -o {out} -m {c_size} -w {c_w} -p {c_p} -f {c_f} -F {c_F} -t {c_threads} --seeds --use-acceleration --print-statistics"
#         for size in data_sizes:
#             for w in w_values:
#                 for p in p_values:
#                     for f in f_values:
#                         for F in F_values:
#                             command = base_command.format(
#                                 pfp = pfp_exe,
#                                 conf =project_base_dir + '/config.ini',
#                                 out =project_base_dir + "/out_w{}_p{}_f{}_F{}_s{}".format(w, p, f, F, size),
#                                 c_size = size, c_w = w, c_p = p, c_f = f, c_F = F, c_threads = n_threads)
#                             out = execute_command(command, 1000000)
#                             with open('out_w{}_p{}_f{}_F{}_s{}.log'.format(w,p,f,F,size), 'wb') as log_file:
#                                 log_file.write(out)
#
#     # Run pfp, without modulo
#     if (RUN_WITHOUT_MODULO):
#         base_command = "/usr/bin/time --verbose {pfp} --configure {conf} -o {out} -m {c_size} -w {c_w} --not-use-modulo -f {c_f} -F {c_F} -t {c_threads} --seeds --use-acceleration --print-statistics"
#         for size in data_sizes:
#             for w in w_values:
#                 for f in f_values:
#                     for F in F_values:
#                         command = base_command.format(
#                             pfp = pfp_exe,
#                             conf =project_base_dir + '/config.ini',
#                             out =project_base_dir + "/out_w{}_f{}_F{}_s{}".format(w, f, F, size),
#                             c_size = size, c_w = w, c_f = f, c_F = F, c_threads = n_threads)
#                         out = execute_command(command, 1000000)
#                         with open('out_w{}_f{}_F{}_s{}.log'.format(w,f,F,size), 'wb') as log_file:
#                             log_file.write(out)

if __name__ == '__main__':
    main_debug()
