# Primer Designer 

Primer Designer is a tool for designing primers for the human reference genome.

It utilises [primer3][primer3-url] and [smalt][smalt-url] to generate primers for given region(s), and outputs a PDF report with the best matched primers.

Primers may be designed for single positions, ranges of regions and also breakpoints for fusion events. All primers are marked up with known SNPs for the given reference from dbSNP and also highlights repeat regions.

Key highlights:
- generates up to 5 pairs of designed primers for given position(s) with GC% and TM
- identifies number of mappings (if not unique)
- highlights SNPs from common SNPs database file
- highlights repeat regions of >6 bases
- highlights bases up and downstream of target

***

## Requirements

- python >=3.6
- [primer3][primer3-url]
- [smalt][smalt-url]
- [samtools][samtools-url]
- GRCh37 & GRCh38 human genome reference FASTA files (+ indexes)
- GRCh37 & GRCh38 common SNPs vcf / vcf.gz (+ indexes)
- Docker

n.b. index files for reference genome and SNPs vcf need generating with both `samtools faidx` and `smalt index INDEX REFSEQ-FILE`

- default setting for smalt index: -k 13 -s 13, where k is word length and s is stepsiz
- recommended setting for 11-24 base reads (i.e. primers): `-k 11 -s 2`, the use of wordlen and stepsiz affects memory requirement, speed, sensitivity and accuracy of mapping, check manual for further details
- smalt index outputs 2 files:
  - `.sma`: compressed set of reference sequences
  - `.smi`: hash index
- smalt manual available here: ftp://ftp.sanger.ac.uk/pub/resources/software/smalt/smalt-manual-0.7.4.pdf
- full required file list:
  - {reference_genome}.fasta
  - {reference_genome}.fasta.fai
  - {reference_genome}.fasta.sma
  - {reference_genome}.fasta.smi
  - {snp}.vcf.gz
  - {snp}.vcf.gz.tbi


***

## General Usage

For usage with Docker, please see below.

Paths to reference files may either be provided through environment variables or passing `--config` and adding paths to the `primer_designer.cfg` file, required vars are given in `example_primer_designer.cfg`.

primer3, smalt and samtools are required to be installed and on path.

Example:

```bash
./primer_designer_region.py -c 1 -p 75761161 --grch37  --config {config-file }# outputs a PDF report around chr 9 pos 12345678
./primer_designer_region.py -c 9 -r 12345678 12346678 --grch37 --config {config-file # outputs a PDF report for a range
./primer_designer_region.py --fusion --b1 9:123456789:b:1 --b2 8:12345678:a:-1 --grch37 --config {config-file # outputs a PDF for fusion design
```
-c
  (required for the position and range only) Specifies the chromosome

-p | -r |
  Specifies the type of

    -p - Specifies position of a base around which primers need to be designed
    -r - Specifies range of nucleotides around which primers need to be designed

--fusion
    Specifies breakpoints locations to design primers for fusion genes; requires passing both `--b1` AND `--b2` each in the format:
         `chr:pos:side:strand`, where side is after (a) or before (b) the breakpoint. The position is included in all of the calculations.
    e.g. `--fusion --b1 5:123456:b:1 --b2 5:456789:a:1`

--grch37 | --grch38
  (required) Specifies the reference genome.

-o
  (optional) Output suffix for naming

-t
  (optional) A flag to output the report in .txt format

-f
  (optional) A flag to change the FLANK region used to design primers around. It is the number of bases before and after the indicated position.

--config (optional) Config file containing paths to required reference files, if not given these MUST be provided as environment variables
<br></br>

***
## Environment Variable
1. `REF_37 / REF_38`: path to human reference genome (fasta) file
2. `SNP_37 / SNP_38`: path to snps database
3. `SNP37_VERSION / SNP38_VERSION`: snp population version e.g. 2.0.1
4. `SNP37_DB / SNP38_DB`: snp population database e.g. gnomad
5. `PRIMER_VERSION`: current primer designer release version

## Docker Usage

A dockerfile is provided that allows for building a full working image of primer designer, this will download and build all necessary tools and has no other requirements for installation other than Docker. The image may be built from the primer_designer dir with: `docker build -t {tag_name} .`

Reference files must still be provided as before, either via a config file or environment variables. It is advised to mount the volume containing these when running and pass the paths relative to the file location within the mounted volume, i.e:

```
docker build -t primer_designer .

docker run -v /home/$USER/reference_files:/reference_files -v $PWD:/home/primer_designer/output --env-file {config-file) {image-name} primer_designer python bin/primer_designer_region.py --chr 12 --pos 56489061 --grch37

docker run -v /home/jason/github/primer_designer/test/reference_files:/reference_files -v /home/jason/github/primer_designer/test/test_output:/home/primer_designer/output --env REF_37=/reference_files/grch37/hs37d5.fa --env SNP_37=/reference_files/grch37/gnomad.genomes.r2.0.1.sites.noVEP.AF-0.01.infoRemoved.vcf.gz --env PRIMER_VERSION=2.0.1 --env SNP37_VERSION=2.0.1 --env SNP37_DB=gnomad primer_designer:2.0.1 python -u bin/primer_designer_region.py --chr 12 --pos 56489061 --grch37

```
In the above example a local dir `/home/$USER/reference_files/` contains the reference files and is mounted in the container at `/reference_files`. The environment variables in `{config-file}` are passed with paths to the files relative from the reference files dir in the container (i.e. `REF37=/reference_files/grch37.fa`).

In addition, the output PDF report will be written to `/home/primer_designer/output/` within the container, and therefore to access this outside the container a volume should be mounted to where you want the report on the local system. To output the report to the current dir it should be passed as in the above example: `-v $PWD:/home/primer_designer/output`
***


## Description
For a detailed description of how to run the script and what has been done for fusion genes please read this [Confluence page][fusion-page-url]

## Authors And Acknowledgement

* Kim Brugger

## Maintainers

Current maintainers:

* ~~Nikita Povarnitsyn~~
* Jethro Rainford (eastgenomics)

## Contributing
Contributions are more than welcome.

[primer3-url]: https://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plusHelp.cgi
[smalt-url]: https://www.sanger.ac.uk/tool/smalt-0/
[fusion-page-url]: https://cuhbioinformatics.atlassian.net/wiki/spaces/BT/pages/481099798/Running+PrimerDesigner+for+fusion+genes
[samtools-url]: http://www.htslib.org/
[dbsnp-url]: https://ftp.ncbi.nih.gov/snp/organisms/
