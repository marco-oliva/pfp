# Prefix Free ParserVCF #
![Build Action Status](https://github.com/marco.oliva/pfp/actions/workflows/cmake.yml/badge.svg)

Tool to build the parse and the dictionary for VCF files using the approach described in Prefix-Free Parsing for Building Big BWTs by Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini.
It produces the same result as running `bigbwt` on the fasta file generated as follow:
```
cat reference.fa | bcftools consensus calls.vcf.gz -H 1 > consensus.fa
```
Symbolic alleles are currently not supported, e.g. `<CN1>`. 

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
Usage: ./pfp++ [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -v,--vcf TEXT ...           List of vcf files. Assuming in genome order!
  -r,--ref TEXT ...           List of reference files. Assuming in genome order!
  -f,--fasta TEXT:FILE        Fasta file to parse.
  -t,--text TEXT:FILE         Text file to parse.
  -o,--out-prefix TEXT        Output prefix
  -m,--max UINT               Max number of samples to analyze
  -w,--window-size UINT:INT in [3 - 200]
                              Sliding window size
  -p,--modulo UINT:INT in [5 - 20000]
                              Module used during parisng
  -j,--threads UINT           Number of threads
  --tmp-dir TEXT:DIR          Tmp file directory
  -c,--compression            Compress the dictionary
  --use-acceleration          Use reference parse to avoid re-parsing
  --print-statistics          Print out csv containing stats
  --version                   Version
  --configure                 Read an ini file
```
