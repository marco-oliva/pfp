# Prefix Free Parser #

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
Usage: pfp++ [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -v,--vcf TEXT ... REQUIRED  List of vcf files. Assuming in genome order
  -r,--ref TEXT ... REQUIRED  List of reference files. Assuming in genome order
  -o,--out-prefix TEXT REQUIRED
                              Output prefix
  -m,--max UINT               Max number of samples to analyze
  -w,--window-size UINT:INT in [3 - 30]
                              Sliding window size
  -p,--module UINT:INT in [50 - 300]
                              Module used during parisng
  -f,--min-frequency FLOAT:FLOAT in [0 - 1]
                              Min frequency for variations
  -t,--threads UINT           Number of threads
  -s,--seeds                  Compute seeded trigger strings
  -c,--compression            Compress the dictionary
  --only-trigger-strings      Generate Only Trigger Strings
  --use-acceleration          Use reference parse to avoid re-parsing
  --print-statistics          Print out csv containing stats
  --version                   Version
  --configure                 Read an ini file
```
