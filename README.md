# Prymer
Command line python programme to aid primer design around genomic fusion breakpoints.

# Details
The primary use of Prymer is to takes a pair of disperate genomic coordinates that represent the breakpoints of a
genomic fusion or structrural variant (SV) and retrieve the genomic sequence immediately up or downstream of each 
breakpoint from the UCSC genome database. The retrieved sequences are then joined and primers are designed to flank 
the SV junction using the primer3-py wrapper for primer3. The top 5 primers from primer 3 are returned.

The output file is a csv file containing the details of the top 5 primer pairs found by primer3. All primers should be checked before use.

Default Primer3 design options are used, but can be controlled using `--p3_sequnce_tags` and/or `--p3_global_tags`, as specified
by the [Primer3 v2.6.1. manual.](https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm)

A single coordinate can also be provided, allowing primers to be designed around a single genomic locus. 


The amount of sequence returned from UCSC can be controlled with the `-l` flag, 

# Example Usage
prymer.py --start_coordinate "chr1:23,456" --end_coordinate "chr6:54,321" --fusion_breakpoint --output_path *path/to/output*

# Options

  -s, --start_coordinate:
    Required. Genomic coordinate for primer design in the format chr1:23456. 
    Commas are allowed, but the full coordinate should be contained by quote marks in this case e.g.
    "chr1:23,456". If --fusion_breakpoint is used, this should be the 5' breakpoint

  -e, --end_coordinate:
    Optional second genomic coordinate for primer design. If *--fusion_breakpoint* is used, this is required and should be the
    3' breakpoint.

  -l, --template_sequence_length:
    Integer length of template sequence in bp to retrieve for primer design. For breakpoint primers, this
    specifies the length of sequence returned either side of the breakpoint. Default 500

  -r, --reference_genome:
    Reference genome to use. Must be a valid genome that can be acessed with UCSC API. Default hg38.",

  -f, --fusion_breakpoint:
    Flag to specify if primers span a breakpoint.

  -o, --output_path:
    Output filepath. Defaults to current directory.

  -n, --output_name:
    Output file name. Default is chrN, where N is the chromosome number passed to start coordinate.

 --start_primer_position:
    Integer to specify the 5' or 3' position of the primer relative to the breakpoint for the 
   *--start_coordinate*. Useful when the sequence of intrest flanking the breakpoints are both on the 
  \- or + strand. Choose one of 5 or 3. Default 5

  --end_primer_position:
    As *--start_primer_position*, for *--end_coordinate*. Default 3.",

  --reverse_complement:
  Specify if the reverse complement of a sequence is required. Options are "start", "end", or "both". 
  For each option the reverse complement of the sequence retrieved up/downstream of the specified coordinate
  will be used for primer design. Useful in conjunction with *--start/end_primer_position*.

  --p3_global_tags:
    Specify additional primer3 global tags as a space separated list in the form "OPTION1_NAME=<value(s)>" "OPTION2_NAME=<value(s)>". 
    The version of primer3 used is 2.6.1. See https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm
    for details of all options available.

 --p3_sequence_tags:
  Specify additional primer3 sequence tags as a space separated list in the form "OPTION1_NAME=<value(s)>" "OPTION2_NAME=<value(s)>".
  The version of primer3 used is 2.6.1. See https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm
  for details of all options available.



