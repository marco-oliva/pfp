#!/usr/bin/env python3

from utils import *

# Set this dir
data_base_dir = '/blue/boucher/marco.oliva/projects/experiments/pfp'
pre_download_data_dir = '/blue/boucher/marco.oliva/data/1kgp'

version_letter = 'a'

# Assuming installed
# - bcftools
# - htslib
# - gzip
# - bgzip

project_base_dir        = os.path.dirname(os.path.abspath(__file__))
date_string             = datetime.datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
data_dir_bigbwt         = "{}/human_tests/{}/bigbwt".format(data_base_dir, date_string)
data_dir_pfp            = "{}/human_tests/{}/pfp".format(data_base_dir, date_string)
common_data_dir         = "{}/human_tests/{}/data".format(data_base_dir, date_string)
common_tools_dir        = "{}/human_tests/tools".format(data_base_dir)
tmp_fasta_dir           = "{}/human_tests/{}/data/tmp".format(data_base_dir, date_string)

# Chromosomes
chromosomes_list = [i for i in range(1, 23)]

w_value     = 10
p_value     = 100
n_threads   = 32

def get_vcf_files(out_dir):
    rootLogger = logging.getLogger()
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
            rootLogger.info('VCF files have to be in {}'.format(pre_download_data_dir + '/vcf/'))
            exit()

        # Filter and index, if not existing
        rootLogger.info('Filtering out symbolyc allels from VCFs')
        execute_command('bcftools view -v snps,indels -m2 -M2 -Oz -o {outv} {inv}'.format(
            inv=out_dir + '/' + chromosome_file_name,
            outv=out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz'), time_it=True)
        execute_command('bcftools index {}'.format(out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz'), time_it=True)

        vcf_files_list.append(out_dir + '/' + chromosome_file_name[:-3] + '.filtered.bgz')
    return vcf_files_list

def get_reference_files(out_dir):
    rootLogger = logging.getLogger()
    out_fasta_list = list()
    for chromosome_id in [str(c) for c in chromosomes_list]:
        if (os.path.exists(pre_download_data_dir + '/reference/' + chromosome_id + '.fa.gz')):
            rootLogger.info('{} already exists'.format(pre_download_data_dir + '/reference/' + chromosome_id + '.fa.gz'))
            execute_command('cp {} {}'.format(pre_download_data_dir + '/reference/' + chromosome_id + '.fa.gz', out_dir + '/' + chromosome_id + '.fa.gz'))
        else:
            rootLogger.info('Reference files have to be in {}'.format(pre_download_data_dir + '/reference/'))
            exit()
        out_fasta_list.append(out_dir + '/' + chromosome_id + '.fa.gz')
    return out_fasta_list

def main():
    logFormatter = logging.Formatter("%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.DEBUG)

    fileHandler = logging.FileHandler("{}/{}_logfile.log".format('.', date_string))
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    parser = argparse.ArgumentParser(description='Testing stuff.')
    parser.add_argument('-s', dest='samples_file', type=str, help='File containing list of samples', required=True)
    parser.add_argument('-t', dest='threads', type=int, help='Number of threads to be used', required=True)
    parser.add_argument('-d', dest='samples_dir', type=str, help='Folder containing pre extracted files', required=True)
    parser.add_argument('--skip-pscan', dest='skip_pscan', help='Skip pscan, run only pfp', action='store_true')
    args = parser.parse_args()

    # Get executables
    mkdir_p(common_tools_dir)
    pfp_exe = get_pfp(common_tools_dir)
    pscan_exe = get_pscan(common_tools_dir)


    if not os.path.exists(args.samples_file):
        rootLogger.info('{} does not exist'.format(args.samples_file))
        return

    with open(args.samples_file) as f_handler:
        samples = f_handler.readlines()
    samples = [x.strip() for x in samples]

    # ============================================================

    # ------------------------------------------------------------

    mkdir_p(data_dir_bigbwt)
    
    out_fasta_multi_sample = data_dir_bigbwt + '/' + str(len(samples)) + '_samples.fa'
    if os.path.exists(out_fasta_multi_sample):
        rootLogger.info('{} already exists, overwriting it'.format(out_fasta_multi_sample))

    multi_sample_file = open(out_fasta_multi_sample, 'w')
    for sample in samples:
        sample_file_path = args.samples_dir + '/' + sample + '/' + sample + '_ALL.fa'
        with open(sample_file_path) as sample_file:
            multi_sample_file.write(sample_file.read())
    multi_sample_file.close()

    # ------------------------------------------------------------
    if (args.skip_pscan):
        # Run pscan.x
        base_command = "{pscan} -t {c_threads} -f {file} -w {window} -p {modulo}"
        command = base_command.format(pscan=pscan_exe, c_threads=n_threads, file=out_fasta_multi_sample, window=w_value,
                                      modulo=p_value)
        execute_command(command, time_it=True)

    else:
        rootLogger.info('Skipped pscan.x')
    # ============================================================

    # Run pfp++
    mkdir_p(data_dir_pfp)

    vcf_files_list = list()
    for chromosome_id in [str(c) for c in range(1,23)]:
        chromosome_file_name = "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5{}.20130502." \
                               "genotypes.vcf.gz".format(chromosome_id, version_letter)

        if (os.path.exists(pre_download_data_dir + '/vcf/' + chromosome_file_name)):
            vcf_files_list.append(pre_download_data_dir + '/vcf/' + chromosome_file_name)
        else:
            rootLogger.info('VCF files have to be in {}'.format(pre_download_data_dir + '/vcf/'))
            exit()


    ref_files_list = list()
    for chromosome_id in [str(c) for c in chromosomes_list]:
        if (os.path.exists(pre_download_data_dir + '/reference/' + chromosome_id + '.fa.gz')):
            ref_files_list.append(pre_download_data_dir + '/reference/' + chromosome_id + '.fa.gz')
        else:
            rootLogger.info('Reference files have to be in {}'.format(pre_download_data_dir + '/reference/'))
            exit()

    pfp_config_file = create_pfp_config_file(vcf_files_list, ref_files_list, data_dir_pfp)
    base_command = "{pfp} --configure {config_file} -t {c_threads} -m {n_samples} " \
                   "--use-acceleration --print-statistics --occurrences -w {window} -p {modulo} -o {out_dir}"
    command = base_command.format(pfp=pfp_exe, c_threads=n_threads, window=w_value,
                                  modulo=p_value, config_file=pfp_config_file, n_samples=len(samples),
                                  out_dir=data_dir_pfp)
    execute_command(command, time_it=True)

if __name__ == '__main__':
    main()
