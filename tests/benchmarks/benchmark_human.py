#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, wget, errno, datetime, logging, gzip, re
from Bio import SeqIO
from multiprocessing import Pool
import tqdm

# Set this dir
data_base_dir = '/blue/boucher/marco.oliva/projects/experiments/pfp'
pre_download_data_dir = '/blue/boucher/marco.oliva/data/1kgp'

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


# Chromosomes
chromosomes_list = [i for i in range(1, 23)]

w_value   = 10
p_value  = 100
n_threads  = 32

#------------------------------------------------------------
# execute command: return command's stdout if everything OK, None otherwise
def execute_command(command, time_it=False, seconds=1000000):
    try:
        if time_it:
            command = '/usr/bin/time --verbose {}'.format(command)
        rootLogger.info("Executing: {}".format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid, stdout=subprocess.PIPE)
        (output, err) = process.communicate()
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        rootLogger.info("Error executing command line:")
        rootLogger.info("\t"+ command)
        return None
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        rootLogger.info("Command exceeded timeout:")
        rootLogger.info("\t"+ command)
        return None
    if output:
        output = output.decode("utf-8")
        rootLogger.info(output)
    if err:
        err = err.decode("utf-8")
        rootLogger.error(err)
    return output


def get_pscan(work_dir):
    if os.path.exists(work_dir + '/Big-BWT'):
        rootLogger.info('Repository already exists')
    else:
        repository = "https://github.com/alshai/Big-BWT.git"
        get_command = "git clone {} {}".format(repository, work_dir + "/Big-BWT")
        execute_command(get_command)
    build_command = "make -C {}".format(work_dir + "/Big-BWT")
    execute_command(build_command)
    return work_dir + '/Big-BWT/newscan.x'

def get_au_pair(work_dir):
    return "binary_path"

def get_pfp(work_dir):
    if os.path.exists('../../build/pfp++'):
        pfp_realpath = os.path.relpath('../../build/pfp++')
        execute_command('cp {} {}'.format(pfp_realpath, work_dir))
        return work_dir + '/pfp++'
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
        chromosome_file_name = "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5{}.20130502." \
                               "genotypes.vcf.gz".format(chromosome_id, version_letter)

        if (os.path.exists(pre_download_data_dir + '/vcf/' + chromosome_file_name)):
            rootLogger.info('Copying {}'.format(chromosome_file_name))
            execute_command('cp {} {}'.format(pre_download_data_dir + '/vcf/' + chromosome_file_name,
                                              out_dir + '/' + chromosome_file_name)
                            )
        else:
            chromosome_url = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/" \
                             + chromosome_file_name
            rootLogger.info('Downloading {}'.format(chromosome_url))
            wget.download(chromosome_url, out_dir + '/' + chromosome_file_name)

        if os.path.exists(pre_download_data_dir + '/vcf/' +  chromosome_file_name[:-3] + '.bgz'):
            rootLogger.info('{} alrady exists'.format(chromosome_file_name[:-3] + '.bgz'))
            execute_command('cp {} {}'.format(pre_download_data_dir + '/vcf/' +  chromosome_file_name[:-3] + '.bgz'),
                            out_dir + '/' + chromosome_file_name[:-3] + '.bgz')
        else:
            if not os.path.exists(out_dir + '/' +  chromosome_file_name[:-3] + '.bgz'):
                execute_command('gunzip {}'.format(out_dir + '/' + chromosome_file_name), time_it=True)
                execute_command('bgzip {}'.format(out_dir + '/' + chromosome_file_name[:-3]), time_it=True)
                execute_command('mv {} {}'.format(out_dir + '/' + chromosome_file_name[:-3] + '.gz',
                                                  out_dir + '/' + chromosome_file_name[:-3] + '.bgz'))

        # Filter and index, if not existing
        rootLogger.info('Filtering out symbolyc allels from VCFs')
        if not os.path.exists(out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz'):
            execute_command('bcftools view -v snps,indels -m2 -M2 -Oz -o {outv} {inv}'.format(
                inv=out_dir + '/' + chromosome_file_name[:-3] + '.bgz',
                outv=out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz'), time_it=True)
        if not os.path.exists(out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz.csi'):
            execute_command('bcftools index {}'.format(out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz'), time_it=True)

        vcf_files_list.append(out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz')
    return vcf_files_list

def get_reference_files(out_dir):
    reference_url = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'
    if (os.path.exists(pre_download_data_dir + '/reference/hs37d5.fa.gz')):
        execute_command('cp {} {}'.format(pre_download_data_dir + '/reference/hs37d5.fa.gz', out_dir + '/hs37d5.fa.gz'))
    else:
        rootLogger.info("Downloading reference files")
        wget.download(reference_url, out_dir + '/' + 'hs37d5.fa.gz')
    rootLogger.info('Extracting chromosomes')

    if not (os.path.exists(out_dir + '/' + 'hs37d5.fa')):
        execute_command("gunzip {}".format(out_dir + '/' + 'hs37d5.fa.gz'))

    out_fasta_list = list()
    for record in SeqIO.parse(out_dir + '/' + 'hs37d5.fa', 'fasta'):
        if (os.path.exists(pre_download_data_dir + '/reference/' + record.id + '.fa.gz')):
            rootLogger.info('{} already exists'.format(pre_download_data_dir + '/reference/' + record.id + '.fa.gz'))
            execute_command('cp {} {}'.format(pre_download_data_dir + '/reference/' + record.id + '.fa.gz', out_dir + '/' + record.id + '.fa.gz'))
        else:
            with gzip.open(out_dir + '/' + record.id + '.fa.gz', 'wt') as output_fasta:
                SeqIO.write(record, output_fasta, "fasta")
        out_fasta_list.append(out_dir + '/' + record.id + '.fa.gz')
    return out_fasta_list

def extract_fasta(out_file_path, ref_file, vcf_file, sample_id):
    bcf_command = "/usr/bin/time --verbose bcftools consensus -H 1 -f {} -s {} {} -o {}".format(
        ref_file, sample_id, vcf_file, out_file_path)
    execute_command(bcf_command, time_it=True)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python ≥ 2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise # nop

def create_pfp_config_file(vcf_list, ref_list, out_dir):
    # Generate config file
    with open(out_dir + '/config.ini', 'w') as config_file:
        config_file.write('# Generated config file\n')
        config_file.write('vcf = {}\n'.format(vcf_list))
        config_file.write('ref = {}\n'.format(ref_list))
    return out_dir + '/config.ini'


def per_sample_fasta(samples_dir, ref_dir, vcf_files_list, sample):
    rootLogger.info('Running on sample: {}'.format(sample))

    sample_dir = samples_dir + '/' + sample
    mkdir_p(sample_dir)
    out_fasta_all_path = sample_dir + '/' + sample + '_ALL.fa'
    if os.path.exists(out_fasta_all_path):
        rootLogger.info('{} already exists, overwriting it'.format(out_fasta_all_path))
    #else:
    out_fasta_all = open(out_fasta_all_path, 'w')
    for idx,vcf_file in enumerate(vcf_files_list):
        chromosome = str(idx + 1)
        ref_fasta = ref_dir + '/' + chromosome + '.fa.gz'
        out_fasta_single_chromosome = sample_dir + '/' + sample + '_' + chromosome + '.fa'
        extract_fasta(out_fasta_single_chromosome, ref_fasta, vcf_file, sample)
        with open(out_fasta_single_chromosome, 'r') as chr_file:
            out_fasta_all.write(chr_file.read())
    out_fasta_all.close()

def main():
    global rootLogger
    logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.DEBUG)

    fileHandler = logging.FileHandler("{0}/{1}.log".format('.', 'logfile'))
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    parser = argparse.ArgumentParser(description='Testing stuff.')
    parser.add_argument('-s', dest='samples_file', type=str, help='File containing list of samples', required=True)
    parser.add_argument('-t', dest='threads', type=int, help='Number of threads to be used', required=True)
    args = parser.parse_args()

    # Get executables
    mkdir_p(common_tools_dir)
    pfp_exe = get_pfp(common_tools_dir)
    pscan_exe = get_pscan(common_tools_dir)

    # Start extracting samples
    if not os.path.exists(args.samples_file):
        rootLogger.info('{} does not exist'.format(args.samples_file))
        return

    with open(args.samples_file) as f_handler:
        samples = f_handler.readlines()
    samples = [x.strip() for x in samples]

    vcf_dir = common_data_dir + '/vcf'
    mkdir_p(vcf_dir)
    vcf_files_list = get_vcf_files(vcf_dir)

    ref_dir = common_data_dir + '/ref'
    mkdir_p(ref_dir)
    ref_files_list = get_reference_files(ref_dir)
    ref_files_list = ref_files_list[0:22]

    samples_dir = common_data_dir + '/samples'
    mkdir_p(samples_dir)

    # ============================================================

    # ------------------------------------------------------------
    # Parallel section
    pool = Pool(processes=args.threads)
    fixed_args = [samples_dir, ref_dir, vcf_files_list]
    execs = []
    for sample in samples:
        execs.append(fixed_args + [sample])

    for _ in tqdm.tqdm(pool.starmap(per_sample_fasta, execs), total=len(execs)):
        pass
    # ------------------------------------------------------------

    out_fasta_multi_sample = samples_dir + '/' + str(len(samples)) + '_samples.fa'
    if os.path.exists(out_fasta_multi_sample):
        rootLogger.info('{} already exists, overwriting it'.format(out_fasta_multi_sample))

    multi_sample_file = open(out_fasta_multi_sample, 'w')
    for sample in samples:
        sample_file_path = samples_dir + '/' + sample + '/' + sample + '_ALL.fa'
        with open(sample_file_path) as sample_file:
            multi_sample_file.write(sample_file.read())
    multi_sample_file.close()

    # ------------------------------------------------------------

    # Run pscan.x
    base_command = "{pscan} -t {c_threads} -f {file} -w {window} -p {modulo}"
    command = base_command.format(pscan=pscan_exe, c_threads=n_threads, file=out_fasta_multi_sample, window=w_value,
                                  modulo=p_value)
    execute_command(command, time_it=True)

    # ============================================================

    # Run pfp++
    mkdir_p(data_dir_pfp)
    pfp_config_file = create_pfp_config_file(vcf_files_list, ref_files_list, data_dir_pfp)
    base_command = "{pfp} --configure {config_file} -t {c_threads} -m {n_samples} " \
                   "--use-acceleration --print-statistics --occurrences -w {window} -p {modulo} -o {out_dir}"
    command = base_command.format(pfp=pfp_exe, c_threads=n_threads, window=w_value,
                                  modulo=p_value, config_file=pfp_config_file, n_samples=len(samples),
                                  out_dir=data_dir_pfp)
    execute_command(command, time_it=True)

if __name__ == '__main__':
    main()
