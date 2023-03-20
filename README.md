# prymer
A python script to automate PCR primer design. Under development

# Usage
prymer.py -c _genomic coordinate(s)_

          -c Required. One or two genomic coordinates in the format "chr1:2345678, chr8:7654321"

          -l Optional. Desired template sequence length in bp. Default 500

          -r Optional. Reference genome to use. Default hg38

          -o Optional. Output file. Default prymer.txt

# Details
The script will use the provided coordinates to retrieve the surrounding genomic sequence from the UCSC database and 
design primers using primer3 (default settings).

If a single coordinate is provided, this will be the center of the retrieved sequence, allowing primers to be designed 
around a single genomic feature. 

Pairs of coordinates should be given in 5' - 3' order relative to the genome, e.g. to design primers surrounding the 
breakpoints of a genomic rearrangement, the 5' breakpoint
of the rearrangement should be specified first, followed by the 3' breakpoint.

Current functionality is basic, and the output file is in a JSON format containing the details of the top 5 primer pairs
found by primer3. All primers should be checked using, e.g. BLAT

# Requirements
requests

primer3-py