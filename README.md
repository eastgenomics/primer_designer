# Primer Designer 

Primer Designer is a program to design primers for the human reference genome.

***

## Usage 

Example:

```bash
./primer_designer_new_strand.py -c 9 -p 12345678 --grch37 #outputs a PDF report around chr 9 pos 12345678 
./primer_designer_new_strand.py -c 9 -r 12345678 12346678 --grch37 #outputs a PDF report for a range
./primer_designer_new_strand.py -b 9:123456789:b:1_8:12345678:a:-1 --grch37 -t #outputs a PDF and TXT reports for a fusion   
```
-c 
  (required for the position and range only) Specifies the chromosome  

-p | -r | -b
  (required) Specifies the type of 

    -p - Specifies position of a base around which primers need to be designed  
    -r - Specifies range of nucleotides around which primers need to be designed
    -b - Specifies breakpoints locations to design primers for fusion genes; requires the input to be in this format:
         chr1:pos:side:strand_chr2:pos:side:strand, where side is after or before the breakpoint. The position is
         included in all of the calculations. 


--hg19 | grch37 | grch38
  (required) Specifies the reference genome. 

-o 
  (optional) Turns on suffix checker

-t 
  (optional) A flag to output the report in .txt format 

-f 
  (optional) A flag to change the FLANK region used to design primers around. It is the numebr of bases before and after the indicated position. 

*** 

## To do
1. Show in report if primers were not found for that particular run 
2. Prevent the script from choosing the best primer near a poly-region 

##Description 
For a detailed description of the script please read this [Confluence page]()

##Authors And Acknowledgement 

* Kim Brugger 

## Maintainers 

Current maintainers: 

* Nikita Povarnitsyn  

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)