# FlexCore

![Python](https://badges.aleen42.com/src/python.svg) ![conda](https://img.shields.io/badge/%E2%80%8B-conda-%2344A833.svg?style=flat&logo=anaconda&logoColor=44A833)

Extract per-site flexible (default ≥95%) bacterial core genomes from read-mapped whole genome alignments. Optionally calculate core genome SNP counts and distances.

For 100% site inclusion core genomes I highly recommend instead using [snp-sites](https://github.com/sanger-pathogens/snp-sites)
for overall speed and performance, followed by [snp-dists](https://github.com/tseemann/snp-dists).

## Dependencies

FlexCore.py is written in python3 and requires the following python packages:

> - numpy
> - screed

## Usage

First set up a conda environment with the appropriate dependencies:

```console
conda env create -f environment.yml
conda activate env-fc
```

Then run the main script `FlexCore.py` referencing your alignment file:

```shell
python FlexCore.py --alignment <alignment file>
```

## Input

Multi-fasta whole genome alignment derived from mapping to a reference and variant calling such as from [snippy](https://github.com/tseemann/snippy), `snippy-core` and `snippy-clean_full_aln`:

```shell
snippy-core --ref ref.fa snippyoutfiles 
snippy-clean_full_aln core.full.aln > clean.full.aln
```

## Output

File | Description
-----|------------
`Coresites.csv` | the core genome alignment  
`rSNPs95.csv`   | the comma separated non-redundant (exclusive) pairwise distances (optional)

> SNP distances are calculated by the number of comparable sites per pair in the core genome alignment. Though a core genome cutoff for the total core alignment may be ≥95%, in principle any number of sites may be missing between a pair of sequences (they may have far fewer sites without gaps or ambiguous bases than the overall alignment size). As such, calculating SNP distance by dividing by the total alignment length may be innapropriate. The SNP counts, SNP distance (SNP count/length of pairwise alignment) and adjusted SNP counts are output (SNP distance times by overall alignment length)

## Options

Flag &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Short flag | Description | Required | Default val
--------------|------------|-------------|----------|--------------
`--alignment` |  `-a` |  Provide path and filename of alignment | ✅
`--cutoff` |     `-c` |  Per-site percent core (integer). Default=95(%) |  | 95
`--nproc` |      `-p` |  Number of processes |                             | 1
`--no-dists` |        |  Do not calculate SNP distances, output core alignment only | | False
`--keepref` |         |  Retain the reference sequence in the core calculation and SNP distances | | False

The `--keepref` option expects the reference file to be named `>Reference` in line with snippy-core output

## Similar software

Other software (that I'm aware of) for extracting flexible core genomes following reference mapping and variant calling:

- <https://github.com/davideyre/runListCompare>  
- <https://github.com/zheminzhou/EToKi> ('align' or 'phylo' module)
