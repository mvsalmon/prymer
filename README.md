# prymer
A python script to automate PCR primer design. Under development

# Usage
prymer.py -c _[options]_

  -h, --help            show this help message and exit      

  -c START_COORDINATE, --start_coordinate START_COORDINATE                                                                                                                  
                        Genomic coordinate for primer design in the format chr1:23456. If --fusion_breakpoint is used, this should be the 5' breakpoint

  -e END_COORDINATE, --end_coordinate END_COORDINATE                                                                                                                        
                        Optional second genomic coordinate for primer design. Default None.If --fusion_breakpoint is used, this is required and should be the 3' breakpoint.

  -l TEMPLATE_SEQUENCE_LENGTH, --template_sequence_length TEMPLATE_SEQUENCE_LENGTH                                                                                          
                        Length of template sequence in bp to use for primer design. Default 500bp  

  -r REFERENCE_GENOME, --reference_genome REFERENCE_GENOME
                        Reference genome to use. Must be valid genome that can be used with UCSC API. Default hg38

  -f, --fusion_breakpoint
                        Specifies if primers must span a fusion breakpoint.

  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Output filepath. Defaults to current directory

  -n OUTPUT_NAME, --output_name OUTPUT_NAME
                        Output file name. Default is chrN, where N is value passed to start coordinate


# Details
The script will use the provided coordinates to retrieve the surrounding genomic sequence from the UCSC database and 
design primers using primer3 (default settings).

If a single coordinate is provided, this will be the center of the retrieved sequence, allowing primers to be designed 
around a single genomic locus. 

Pairs of coordinates should be given in 5' - 3' order relative to the genome, e.g. to design primers surrounding the 
breakpoints of a genomic rearrangement, the 5' breakpointof the rearrangement should be specified first, followed by 
the 3' breakpoint.

Current functionality is basic, and the output file is a csv file containing the details of the top 5 primer pairs
found by primer3. All primers should be checked before use.

# Requirements
requests~=2.31.0
primer3-py~=1.0.0
pandas~=1.5.3