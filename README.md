# Primer Designer 

Primer Designer is a program to design primers for the human reference genome.

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
For a detailed description of how to run the script and what has been done for fusion genes please read this [Confluence page](https://cuhbioinformatics.atlassian.net/wiki/spaces/BT/pages/481099798/Running+PrimerDesigner+for+fusion+genes)

## Authors And Acknowledgement 

* Kim Brugger 

## Maintainers 

Current maintainers: 

* ~~Nikita Povarnitsyn~~ 
* Jethro Rainford (eastgenomics)

## Contributing
Contributions are more than welcome. 