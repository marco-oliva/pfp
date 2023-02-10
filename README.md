# Prefix-Free Parsing #
[![CMake](https://github.com/marco-oliva/pfp/actions/workflows/cmake.yml/badge.svg?branch=master)](https://github.com/marco-oliva/pfp/actions/workflows/cmake.yml)
[![Conda](https://anaconda.org/bioconda/pfp/badges/version.svg)](https://anaconda.org/bioconda/pfp)

Tool to build the parse and the dictionary for VCF files using the approach described in Prefix-Free Parsing for Building Big BWTs by Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini.
It produces the same result as running `bigbwt` on the fasta file generated as follow:
```
cat reference.fa | bcftools consensus calls.vcf.gz -H 1 > consensus.fa
```
Symbolic alleles are currently not supported, e.g. `<CN1>`. 

### Bioconda ###
PFP is available on `bioconda`:

```bash
conda install -c bioconda -c conda-forge pfp
pfp++ --help
```

### Docker ###
PFP is available on docker:

```bash
docker pull moliva3/pfp:latest
docker run moliva3/pfp:latest pfp++ --help
```

If using singularity:
```bash
singularity pull pfp_sif docker://moliva3/pfp:latest
./pfp_sif pfp++ --help
```

### Build ###

#### Dependencies ####

* Htslib
* OpenMP

#### Build Instructions ####

```
git clone https://github.com/marco-oliva/pfp.git
cd pfp
mkdir build && cd build
cmake ..
make
```

### Usage ###

```
PFP++
Usage: pfp++ [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -v,--vcf TEXT ...           List of comma ',' separated vcf files. Assuming in genome order!
  -r,--ref TEXT ...           List of comma ',' separated reference files. Assuming in genome order!
  -f,--fasta TEXT:FILE        Fasta file to parse.
  -i,--int32t TEXT:FILE       Integers file to parse.
  --int-shift INT:INT in [0 - 200]
                              Each integer i in int32t input is interpreted as (i + int-shift).
  -H,--haplotype TEXT         Haplotype. [1,2,12]
  -t,--text TEXT:FILE         Text file to parse.
  -o,--out-prefix TEXT        Output prefix.
  -m,--max UINT               Max number of samples to analyze.
  -S,--samples TEXT           File containing the list of samples to parse.
  -w,--window-size UINT:INT in [3 - 200]
                              Sliding window size.
  -p,--modulo UINT:INT in [5 - 20000]
                              Modulo used during parisng.
  -j,--threads UINT           Number of threads.
  --tmp-dir TEXT:DIR          Temporary files directory.
  -c,--compression            Also output compressed the dictionary.
  --use-acceleration          Use reference parse to avoid re-parsing.
  --print-statistics          Print out csv containing stats.
  --verbose                   Verbose output.
  --version                   Version number.
  --configure                 Read an ini file.
```
