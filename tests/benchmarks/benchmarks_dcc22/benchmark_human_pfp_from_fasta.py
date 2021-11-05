#!/usr/bin/env python3

from utils import *

# Set this dir
data_base_dir = '/blue/boucher/marco.oliva/projects/experiments/pfp/DCC22/vcf_to_fa'

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
chromosomes_list = [17, 18, 19]

w_value     = 10
p_value     = 100

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
    parser.add_argument('-d', dest='samples_dir', type=str, help='Folder containing pre extracted files', required=True)
    parser.add_argument('-m', default=0, dest='max_samples', type=int, help='Out of the samples list, get first m', required=True)
    args = parser.parse_args()

    # Get executables
    mkdir_p(common_tools_dir)
    pfp_exe = get_pfp(common_tools_dir)

    if not os.path.exists(args.samples_file):
        rootLogger.info('{} does not exist'.format(args.samples_file))
        return

    with open(args.samples_file) as f_handler:
        samples = f_handler.readlines()
    samples = [x.strip() for x in samples]
    if (args.max_samples != 0):
        samples = samples[0 : args.max_samples]

    # ============================================================

    mkdir_p(data_dir_bigbwt)

    # ------------------------------------------------------------
    out_fasta_multi_sample = data_dir_bigbwt + '/' + str(len(samples)) + '_samples.fa'
    if os.path.exists(out_fasta_multi_sample):
        rootLogger.info('{} already exists, overwriting it'.format(out_fasta_multi_sample))

    multi_sample_file = open(out_fasta_multi_sample, 'w')
    for idx, sample in enumerate(samples):
        sample_file_path = args.samples_dir + '/' + sample + '/' + sample + '_ALL_H1_H2.fa'
        rootLogger.info('Cpying [{}/{}] {}'.format(idx, args.max_samples,sample_file_path))
        with open(sample_file_path) as sample_file:
            multi_sample_file.write(sample_file.read())
    multi_sample_file.close()

    # Run pfp++ from fasta
    base_command = "{pfp} -w {window} -p {modulo} -f {input_file}"
    command = base_command.format(pfp=pfp_exe, window=w_value, modulo=p_value, input_file=out_fasta_multi_sample)
    execute_command(command, time_it=True)

    # ============================================================

if __name__ == '__main__':
    main()
