# TODO allow coordinate input with commas eg 1,234,567
# TODO take file with coordinates as input?
# TODO add PRIMER_EXPLAIN info to output file

import argparse
import re

import requests
import pandas as pd
import primer3
from primer3 import p3helpers

class Primer:
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
        self.start_primer_position = str(args.start_primer_position)
        self.end_primer_position = str(args.end_primer_position)
        self.rev_comp = args.reverse_complement

        # variables to store primer3 output
        self.primer3_info = {}
        self.primer3_pairs = {}
        self.pair_explain = ""

        # primer3 options
        self.seq_args = {"SEQUENCE_ID": self.output_name, "SEQUENCE_TEMPLATE": None}
        if not args.primer3_global_opts:
            self.global_args = {}
        else:
            self.global_args = self._parse_primer3_global_opts(args.primer3_global_opts)

        #TODO Validate output path
        #TODO check for existing file?

        self._run()

    def _parse_primer3_global_opts(self, p3_global_opts):
        global_opts_dict = {}
        for option in p3_global_opts:
            param, value = option.split(":")
            global_opts_dict[param.strip()] = value.strip()

        return global_opts_dict

    def _run(self):
        """Main program control"""
        # single or pair of non-fusion (i.e. on same chrom) coordinate primers
        # TODO set include region for two non-fusion primers
        if not self.fusion_breakpoint:
            # raise exception if both coordinates are on different chroms.
            try:
                self.UCSC_response = self.UCSC_request(
                    self.start_coordinate, self.end_coordinate
                )
            except ChromMismatch as error:
                print(error)
                exit(1)
            # store template sequence from API request
            self.template_sequence = self.UCSC_response["dna"]
            self.seq_args["SEQUENCE_TEMPLATE"] = self.template_sequence
            # print(self.UCSC_response["dna"])
            # design primers
            self.primers = self.design_primers()

            self.write_output()

        # fusion breakpoint primers. Must be a pair. Call UCSC API request for each breakpoint.
        if self.fusion_breakpoint and self.end_coordinate:
            self.UCSC_start_breakpoint_response = self.UCSC_request(
                self.start_coordinate, breakpoint_position=self.start_primer_position
            )
            # print(self.UCSC_start_breakpoint_response)
            self.UCSC_end_breakpoint_response = self.UCSC_request(
                self.end_coordinate, breakpoint_position=self.end_primer_position
            )
            # print(self.UCSC_end_breakpoint_response)

            # concatenate sequence at each breakpoint
            self.breakpoint_sequence = self._build_breakpoint()

            # update primer3 options
            self.seq_args["SEQUENCE_TEMPLATE"] = self.breakpoint_sequence
            # include breakpoint region, 100bp either side
            self.seq_args["SEQUENCE_TARGET"] = [self.seq_len - 100, 200]
            self.primers = self.design_primers()

            # write output
            self.write_output()

    def UCSC_request(
        self, start_coordinate, end_coordinate=None, breakpoint_position=None
    ):
        """Handles different sets of coordinates (eg, single, pair, fusion), then
        parses and requests sequence data from UCSC. Returns json."""
        print("INFO: Querying UCSC database...")
        # two coordinates on same chromosome
        if end_coordinate:
            end_chrom, end = end_coordinate.split(":")
            start_chrom, start = start_coordinate.split(":")
            # check both coords are on the same chromosome and fusion primers are not required
            if start_chrom != end_chrom:
                raise ChromMismatch(
                    f"Error! {start_coordinate} and {end_coordinate} are on different chromosomes. "
                    f"Specify the same chromosome or use --fusion_breakpoint."
                )
            url = f"https://api.genome.ucsc.edu/getData/sequence?genome={self.ref_genome};chrom={start_chrom};start={start};end={end}"
        # single coordinate
        elif not breakpoint_position:
            # request for single coordinate
            chrom, start = start_coordinate.split(":")
            # given coordinate will be in the center of the returned sequence,
            # so adjust start and end position by seq_len/2
            flanking_len = round(self.seq_len / 2)
            end = int(start) + flanking_len
            start = int(start) - flanking_len
            url = f"https://api.genome.ucsc.edu/getData/sequence?genome={self.ref_genome};chrom={chrom};start={start};end={end}"
        # fusion breakpoints
        else:
            # adjust start or end position for 5'/3' break as required
            if breakpoint_position == "5":
                chrom, end = start_coordinate.split(":")
                start = int(end) - self.seq_len
            elif breakpoint_position == "3":
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

    def _build_breakpoint(self):
        """concantenate breakpoint sequences"""
        # reverse complement sequences as required
        if self.rev_comp == "start":
            self.UCSC_start_breakpoint_response["dna"] = self._reverse_comp(
                self.UCSC_start_breakpoint_response["dna"]
            )
        if self.rev_comp == "end":
            self.UCSC_end_breakpoint_response["dna"] = self._reverse_comp(
                self.UCSC_end_breakpoint_response["dna"]
            )
        breakpoint_sequence = "".join(
            [
                self.UCSC_start_breakpoint_response["dna"],
                self.UCSC_end_breakpoint_response["dna"],
            ]
        )
        return breakpoint_sequence

    # p3helpers only included in v1.2.0+
    def _reverse_comp(self, sequence):
       return primer3.p3helpers.reverse_complement(sequence)

    def design_primers(self):
        """design PCR primers using primer3 with default options"""
        # TODO primer design options?
        print("INFO: Designing primers...")
        primers = primer3.design_primers(
            seq_args=self.seq_args, global_args=self.global_args
        )
        parsed_primers = self._parse_primer3(primers)
        self.pair_explain = primers['PRIMER_PAIR_EXPLAIN']
        return parsed_primers

    def _design_breakpoint_primers(self):
        """placeholder for possible future fusion stuff"""
        # primers = self._design_primers(self.breakpoint_sequence_template)
        # return primers
        pass

    def _parse_primer3(self, primer3_output):
        """parse primer3 output dictionary to be more manageable."""
        print("INFO: Parsing results...")
        parsed = {}
        # store primer3 output as a dict. Keys are pair IDs, values are dict of results for that pair
        # compatible with primer3-py v2.0.0 output by skipping PRIMER_{PAIR, LEFT, RIGHT, INTERNAL} list in output dict
        # should also be backwards compatible with older versions
        for key, value in primer3_output.items():
            # catch index error from PRIMER_{PAIR, LEFT, RIGHT, INTERNAL} list to skip over
            try:
                pair_id = str("PAIR_" + key.split(sep="_")[2])
                if pair_id not in parsed:
                    parsed[pair_id] = {}
                # replace any numbers in key string for ease of output formatting
                key = re.sub(r"_\d", "", key)
                parsed[pair_id][key] = value
            except IndexError:
                pass

        # separate primer3 info and per-pair results
        for pair, result in parsed.items():
            if pair in ["PAIR_EXPLAIN", "PAIR_NUM"]:
                self.primer3_info[pair] = result
            else:
                self.primer3_pairs[pair] = result

    def write_output(self):
        """save results to disk"""
        print("INFO: Writing output...")
        outpath = f"{self.output_path}/{self.output_name}.csv"
        pd.DataFrame.from_dict(self.primer3_pairs).to_csv(outpath)
        # Show primer design details - useful if no pairs returned
        print(f"Design details: {self.pair_explain}")


# exceptions
class ChromMismatch(Exception):
    """raise this if start and end chroms are different in non-breakpoint situation"""

    pass


def prymer_main():
    parser = argparse.ArgumentParser(
        description="Design PCR primers for given genomic coordinates"
    )
    parser.add_argument(
        "-s",
        "--start_coordinate",
        help="Genomic coordinate for primer design in the format chr1:23456. "
        "If --fusion_breakpoint is used, this should be the 5' breakpoint",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--end_coordinate",
        help="Optional second genomic coordinate for primer design. Default None."
        "If --fusion_breakpoint is used, this is required and should be the "
        "3' breakpoint.",
        default=None,
    )
    parser.add_argument(
        "-l",
        "--template_sequence_length",
        help="Length of template sequence in bp to use for primer design. For breakpoint primers, this"
        "specifies the length of sequence returned either side of the breakpoint. Default 500bp",
        type=int,
        default=500,
    )
    parser.add_argument(
        "-r",
        "--reference_genome",
        help="Reference genome to use. Must be valid genome that can be used with UCSC API. Default hg38.",
        default="hg38",
    )
    parser.add_argument(
        "-f",
        "--fusion_breakpoint",
        help="Specifies if primers must span a fusion breakpoint.",
        action="store_true",
    )
    # TODO change to path object/verify path is valid
    parser.add_argument(
        "-o",
        "--output_path",
        help="Output filepath. Defaults to current directory",
        default="./",
    )
    parser.add_argument(
        "-n",
        "--output_name",
        help="Output file name. Default is chrN, where N is value passed to start coordinate",
        default=None,
    )
    parser.add_argument(
        "--start_primer_position",
        help="""Specify the relative ("5'" or "3'") position of the primer for the 
                        --start_coordinate. Should be used when a each side of a breakpoint are both on the 
                        - or + strand. Choose one of 5 or 3. Default 5. Note: do not include the trailing "'".""",
        default="5",
    )
    parser.add_argument(
        "--end_primer_position",
        help="As --start_primer_position, for --end_coordinate. Default 3.",
        default="3",
    )
    parser.add_argument(
        "--reverse_complement",
        help="Specify if either the reverse complement of either the start or end sequence is required."
        'Options are "start" or "end"',
        choices=["start", "end"],
    )
    parser.add_argument(
        "--primer3_global_opts",
        nargs="*",
        help="""Specify additional primer3 global options in the form OPTION_NAME:<value>. The version of primer3 used 
                is 2.6.1. See https://htmlpreview.github.io/?https://github.com/primer3-org/primer3/blob/v2.6.1/src/primer3_manual.htm
                for details of all options available."""
    )

    args = parser.parse_args()
    if args.output_name is None:
        args.output_name = args.start_coordinate.split(sep=":")[0]

    # print("Designing primers...")
    Primer(args)

    # return primers


if __name__ == "__main__":
    prymer_main()
