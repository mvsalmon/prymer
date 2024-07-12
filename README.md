# Prymer
Command line python programme to aid primer design around genomic structural variant (SV) breakpoints.

# Details
The primary use of Prymer is to takes a pair of genomic coordinates that represent the breakpoints of an SV and 
retrieve the genomic sequence immediately up of downstream of the breakpoints from the UCSC genome database. The
retrieved sequences are then joined and primers are designed (using the primer3-py wrapper for primer3) to flank the
SV junction, returning the top 5 primers from primer 3.

The amount of sequence returned from UCSC can be controlled with the `-l` flag, and primer 3 options can be specified
using `--p3_sequnce_tags` and/or `--p3_global_tags`.

A single coordinate can also be provided, allowing primers to be designed around a single genomic locus. 






Current functionality is basic, and the output file is a csv file containing the details of the top 5 primer pairs
found by primer3. All primers should be checked before use.

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


