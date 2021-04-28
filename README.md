# Primer Designer 

Primer Designer is a tool for designing primers for the human reference genome.

It utilises [primer3][primer3-url] and [smalt][smalt-url] to generate primers for given region(s), and outputs a PDF report with the best matched primers.

Primers may be designed for single positions, ranges of regions and also breakpoints for fusion events. All primers are marked up with known SNPs for the given reference from dbSNP and also highlights repeat regions.

***

## Requirements

- python >=3.6 
- [primer3][primer3-url] (2.3.7)
- [smalt][smalt-url] (0.7.6)
- [samtools][samtools-url] (1.5)
- GRCh37 & GRCh38 reference FASTA files
- GRCh37 & GRCh38 [dbSNP VCF][dbsnp-url]
- Docker (optional)

***

## General Usage 

For useage with Docker, please see below.

Paths to reference files may either be provided through environment variables or passing `--config` and adding paths to the `primer_designer.cfg` file, required vars are given in `example_primer_designer.cfg`.

primer3, smalt and samtools are required to be installed and on path.

Example:

```bash
./primer_designer_region.py -c 1 -p 75761161 --grch37 #outputs a PDF report around chr 9 pos 12345678 
./primer_designer_region.py -c 9 -r 12345678 12346678 --grch37 #outputs a PDF report for a range
./primer_designer_region.py -b 9:123456789:b:1_8:12345678:a:-1 --grch37 -t #outputs a PDF and TXT reports for a fusion   
```
-c 
  (required for the position and range only) Specifies the chromosome  

-p | -r | -b
  (required) Specifies the type of 

    -p - Specifies position of a base around which primers need to be designed  
    -r - Specifies range of nucleotides around which primers need to be designed
    -b - Specifies breakpoints locations to design primers for fusion genes; requires the input to be in this format:
         chr1:pos:side:strand_chr2:pos:side:strand, where side is after (a) or before (b) the breakpoint. The position is included in all of the calculations. 

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

Primer designer may also be called via the `bulk_design.py` script, this allows for passing a `.txt` file of multiple regions at once and outputs a `.zip` file of designer primer reports.
*** 

## Docker Usage

A dockerfile is provided that allows for building a full working image of primer designer, this will download and build all necessary tools and has no other requirements for installation other than Docker. The image may be built from the primer_designer dir with: `docker build -t {tag_name} .`

Reference files must still be provided as before, either via a config file or environment variables. It is advised to mount the volume containing these when running and pass the paths relative to the file location within the mounted volume, i.e:

```
# running primer designer from docker
docker run 
  -v /home/$USER/reference_files:/reference_files 
  -v $PWD:/home/primer_designer/output 
  --env REF_37=/reference_files/grch37/human_g1k_v37.fasta 
  --env DBSNP_37=/reference_files/dbsnp_37.vcf.gz 
  primer_designer 
  python3 primer_designer/primer_designer_region.py -c 8 -p 21988118 --grch37
```
In the above example a local dir `reference_files/` contains the files and is mounted in the container at `/reference_files`. The environment variables `REF_37` and `DBSNP_37` are passed with paths to the files relative from the reference files dir.

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
[db-snp]: https://ftp.ncbi.nih.gov/snp/organisms/
