# TODO update primer output parsing for compatability with primer3-py v2.0.0

import argparse
import re

import requests
import primer3
import pandas as pd

class Primer():
    """Class to manage primer objects.
    seq_len = length of sequence in bp that will be returned from UCSC query. Value will be used to calculate the start
    and end positions of the sequence requested from UCSC.
    When a single coordinate is supplied, this will be at the center of the returned sequence.
    For two coordinates (i.e. fusion), each coordinate should represent the breakpoint, and will be either the
    start or end position of the returned sequence."""

    def __init__(self, args):
        self.start_coordinate = args.start_coordinate
        self.end_coordinate = args.end_coordinate
        self.ref_genome = args.reference_genome
        self.seq_len = args.template_sequence_length
        self.fusion_breakpoint = args.fusion_breakpoint
        self.output_path = args.output_path
        self.output_name = args.output_name

        # variables to store primer3 output
        self.primer3_info = {}
        self.primer3_pairs = {}

        # primer3 options
        self.seq_args = {'SEQUENCE_ID': self.output_name, 'SEQUENCE_TEMPLATE': None}
        self.global_args = {}

        self._run()


    def _run(self):
        """Main program control"""
        # single or pair of non-fusion (i.e. on same chrom) coordinate primers
        if not self.fusion_breakpoint:
            try:
                # raise exception if both coordinates are on different chroms.
                self.UCSC_response = self.UCSC_request(self.start_coordinate, self.end_coordinate)
            except ChromMismatch as error:
                print(error)
                exit(1)
            # store template sequence from API request
            #self.template_sequence = self.UCSC_response["dna"]
            self.seq_args['SEQUENCE_TEMPLATE'] = self.template_sequence
            # print(self.UCSC_response["dna"])
            # design primers
            self.primers = self.design_primers()


        # fusion breakpoint primers. Must be a pair. Call UCSC API request for each breakpoint.
        if self.fusion_breakpoint and self.end_coordinate:
            self.UCSC_start_breakpoint_response = self.UCSC_request(self.start_coordinate, breakpoint_position="5'")
            #print(self.UCSC_start_breakpoint_response)
            self.UCSC_end_breakpoint_response = self.UCSC_request(self.end_coordinate, breakpoint_position="3'")
            #print(self.UCSC_end_breakpoint_response)

            # concatenate sequence at each breakpoint
            self.breakpoint_sequence = self._build_breakpoint()

            # update primer3 options
            self.seq_args['SEQUENCE_TEMPLATE'] = self.breakpoint_sequence
            # include breakpoint region, 100bp either side
            self.seq_args['SEQUENCE_TARGET'] = [self.seq_len - 100, 200]
            self.primers = self.design_primers()

            #write output
            outpath = f'{self.output_path}{self.output_name}.csv'
            pd.DataFrame.from_dict(self.primer3_pairs).to_csv(outpath)
        # print(pd.DataFrame.from_dict(self.primers))

    def UCSC_request(self, start_coordinate, end_coordinate=None, breakpoint_position=None):
        """Handles different sets of coordinates (eg, single, pair, fusion), then
         parses and requests sequence data from UCSC. Returns json."""
        # two coordinates on same chromosome
        if end_coordinate:
            end_chrom, end = end_coordinate.split(":")
            start_chrom, start = start_coordinate.split(":")
            # check both coords are on the same chromosome and fusion primers are not required
            if start_chrom != end_chrom:
                raise ChromMismatch(f"Error! {start_coordinate} and {end_coordinate} are on different chromosomes. "
                                    f"Specify the same chromosome or use --fusion_breakpoint.")

            url = f"https://api.genome.ucsc.edu/getData/sequence?genome={self.ref_genome};chrom={start_chrom};start={start};end={end}"

        elif not breakpoint_position:
            # request for single coordinate
            chrom, start = start_coordinate.split(":")
            # given coordinate will be in the center of the returned sequence,
            # so adjust start and end position by seq_len/2
            flanking_len = round(self.seq_len/2)
            end = int(start) + flanking_len
            start = int(start) - flanking_len
            url = f"https://api.genome.ucsc.edu/getData/sequence?genome={self.ref_genome};chrom={chrom};start={start};end={end}"

        else:
            # request for breakpoint primers
            # adjust start or end position for 5'/3' break as required
            if breakpoint_position == "5'":
                chrom, end = start_coordinate.split(":")
                start = int(end) - self.seq_len
            elif breakpoint_position == "3'":
                chrom, start = start_coordinate.split(":")
                end = int(start) + self.seq_len

            url = f"https://api.genome.ucsc.edu/getData/sequence?genome={self.ref_genome};chrom={chrom};start={start};end={end}"

        # make API request and check ok
        response = requests.get(url)
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as error:
            print(error)
            exit(1)

        # print(response.url)
        return response.json()

    def design_primers(self):
        """design PCR primers using primer3 with default options"""
        # TODO primer design options?
        # TODO set include breakpoint region for fusion primer
        primers = primer3.design_primers(
            seq_args=self.seq_args,
            global_args=self.global_args)
        parsed_primers = self._parse_primer3(primers)
        return parsed_primers

    def _build_breakpoint(self):
        '''concantenate breakpoint sequences'''
        breakpoint_sequence = ''.join([self.UCSC_start_breakpoint_response['dna'],
                                       self.UCSC_end_breakpoint_response['dna']])
        return breakpoint_sequence

    def _design_breakpoint_primers(self):
        """placeholder for possible future fusion stuff"""
        #primers = self._design_primers(self.breakpoint_sequence_template)
        #return primers
        pass
    def _parse_primer3(self, primer3_output):
        """parse primer3 output dictionary to be more manageable."""
        parsed = {}
        # store primer3 output as a dict. Keys are pair IDs, values are dict of results for that pair
        for key, value in primer3_output.items():
            pair_id = str('PAIR_' + key.split(sep="_")[2])
            if pair_id not in parsed:
                parsed[pair_id] = {}
            # replace any numbers in key string for ease of output formatting
            key = re.sub('_\d', '', key)
            parsed[pair_id][key] = value

        # separate primer3 info and per-pair results
        for pair, result in parsed.items():
            if pair in ['PAIR_EXPLAIN', 'PAIR_NUM']:
                self.primer3_info[pair] = result
            else:
                self.primer3_pairs[pair] = result

# exceptions
class ChromMismatch(Exception):
    '''raise this if start and end chroms are different in non-breakpoint situation'''
    pass

def prymer_main():
    parser = argparse.ArgumentParser(description='Design PCR primers for given genomic coordinates')
    # TODO add output file name/path options
    parser.add_argument('-s',
                        '--start_coordinate',
                        help="Genomic coordinate for primer design in the format chr1:23456. "
                             "If --fusion_breakpoint is used, this should be the 5' breakpoint",
                        required=True)
    parser.add_argument('-e',
                        '--end_coordinate',
                        help="Optional second genomic coordinate for primer design. Default None."
                             "If --fusion_breakpoint is used, this is required and should be the " 
                             "3' breakpoint.",
                        default=None)
    parser.add_argument('-l',
                        '--template_sequence_length',
                        help='Length of template sequence in bp to use for primer design. For breakpoint primers, this'
                             'specifies the length of sequence returned either side of the breakpoint. Default 500bp',
                        type=int,
                        default=500)
    parser.add_argument('-r',
                        '--reference_genome',
                        help='Reference genome to use. Must be valid genome that can be used with UCSC API. Default hg38',
                        default='hg38')
    parser.add_argument('-f',
                        '--fusion_breakpoint',
                        help="Specifies if primers must span a fusion breakpoint.",
                        action="store_true")
    parser.add_argument('-o',
                        '--output_path',
                        help="Output filepath. Defaults to current directory",
                        default="./")
    parser.add_argument('-n',
                        '--output_name',
                        help="Output file name. Default is chrN, where N is value passed to start coordinate",
                        default=None)

    args = parser.parse_args()
    if args.output_name is None:
        args.output_name = args.start_coordinate.split(sep=":")[0]

    print("Designing primers...")
    Primer(args)

    # return primers

if __name__ == '__main__':
    prymer_main()
