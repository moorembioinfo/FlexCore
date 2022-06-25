# FlexCore
Extract per-site flexible (default ≥95%) bacterial core genomes from read-mapped whole genome alignments. Optionally calculate core genome SNP counts and distances.

For 100% site inclusion core genomes I highly recommend instead using snp-sites (https://github.com/sanger-pathogens/snp-sites)
for overall speed and performance, followed by snp-dists (https://github.com/tseemann/snp-dists). 

## Dependencies

FlexCore.py is written in python3 and requires the following python packages: 
>os; sys; argparse; random; numpy; pandas; itertools; concurrent

## Usage
	python FlexCore.py --alignment [alignment file]

## Input
Multi-fasta whole genome alignment derived from mapping to a reference and variant calling such as from snippy, snippy-core and snippy-clean_full_aln:

	snippy-core --ref ref.fa snippyoutfiles 
	snippy-clean_full_aln core.full.aln > clean.full.aln
https://github.com/tseemann/snippy

## Output
*Coresites.csv*, the core genome alignment  
*rSNPs95.csv*, the comma separated non-redundant (exclusive) pairwise distances (optional)

>SNP distances are calculated by the number of comparable sites per pair in the core genome alignment. Though a core genome cutoff for the total core alignment may be ≥95%, in principle any number of sites may be missing between a pair of sequences (they may have far fewer sites without gaps or ambiguous bases than the overall alignment size). As such, calculating SNP distance by dividing by the total alignment length may be innapropriate. The SNP counts, SNP distance and adjusted SNP counts are output (SNP distance times by overall alignment length)


## Options


	'--alignment', '-a', help='Provide path and filename of alignment', required=True
	'--cutoff', '-c', help='Per-site percent core (integer). Default=95(%)', type=int, default=95
	'--nproc', '-p', help='Number of processes; type=int, default=1
    '--nodists', help='Dont calculate SNP distances, only output core alignment; action='store_true'
	'--keepref', help='Retain the reference sequence in the core calculation and SNP distances; action='store_true' 

The --keepref option expects the reference file to be named '>Reference' in line with snippy-core output
