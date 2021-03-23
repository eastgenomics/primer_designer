# Primer Designer 

Primer Designer is a tool for designing primers for the human reference genome.

It utilises [primer3][primer3-url] and [smalt][smalt-url] to generate primers for given region(s), and outputs a PDF report with the best matched primers.

Primers may be designed for single positions, ranges of regions and also breakpoints for fusion events. All primers are marked up with known SNPs for the given reference from dbSNP and also highlights repeat regions.

***

## Usage 

`config.py` must be first populated with paths to required tools and files.

Example:

```bash
./primer_designer_region.py -c 9 -p 12345678 --grch37 #outputs a PDF report around chr 9 pos 12345678 
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
  (optional) Turns on suffix checker

-t 
  (optional) A flag to output the report in .txt format 

-f 
  (optional) A flag to change the FLANK region used to design primers around. It is the number of bases before and after the indicated position. 

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
