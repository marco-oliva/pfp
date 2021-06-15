import sys, time, argparse, subprocess, os, wget, errno, datetime, logging, gzip, re
from Bio import SeqIO
from multiprocessing import Pool
import tqdm


#------------------------------------------------------------
# execute command: return command's stdout if everything OK, None otherwise
def execute_command(command, time_it=False, seconds=1000000):
    try:
        if time_it:
            command = '/usr/bin/time --verbose {}'.format(command)
        rootLogger.info("Executing: {}".format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = process.communicate()
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        rootLogger.info("Error executing command line")
        return None
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        rootLogger.info("Command exceeded timeout")
        return None
    if output:
        output = output.decode("utf-8")
        rootLogger.info(output)
    if err:
        err = err.decode("utf-8")
        rootLogger.info("\n" + err)
    return output

#------------------------------------------------------------
# download and compile pscan
def get_pscan(work_dir):
    if os.path.exists(work_dir + '/Big-BWT/newscan.x'):
        rootLogger.info('newscan.x already exists')
    else:
        repository = "https://github.com/alshai/Big-BWT.git"
        get_command = "git clone {} {}".format(repository, work_dir + "/Big-BWT")
        execute_command(get_command)
        build_command = "make -C {}".format(work_dir + "/Big-BWT")
        execute_command(build_command)
    return work_dir + '/Big-BWT/newscan.x'


#------------------------------------------------------------
# download and compile aupair
def get_au_pair(work_dir):
    return "binary_path"

#------------------------------------------------------------
# download and compile pfp++
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

#------------------------------------------------------------
# extract fasta using bcftools
def extract_fasta(out_file_path, ref_file, vcf_file, sample_id):
    bcf_command = "bcftools consensus -H 1 -f {} -s {} {} -o {}".format(
        ref_file, sample_id, vcf_file, out_file_path)
    execute_command(bcf_command, time_it=True)

#------------------------------------------------------------
# mkdir
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


