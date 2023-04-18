# TODO command line argument handling
# TODO handle multiple coordinates for fusions
# TODO error handling

import argparse

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

    def __init__(self, start_pos, end_pos = None, ref_genome = 'hg38', seq_len = 500):
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.ref_genome = ref_genome
        self.seq_len = seq_len

        # self.coordinates = args.coordinates
        # self.ref_genome = args.reference_genome
        # self.template_sequence_length = args.template_sequence_length

        if self.end_pos is None:
            try:
                self.sequence = self._UCSC_request()
            except ChromMismatch as error:
                print(error)


            self.parsed_coordinate = self._parse_coordinate(self.coordinates[0])
            self.sequence_data = self._UCSC_request(self.parsed_coordinate)
            self.primers = self._design_primers(self.sequence_data['dna'])
            print(self.primers)
            print(pd.DataFrame.from_dict(self.primers))

        elif len(self.coordinates) == 2:
            self.parsed_left_coordinate = self._parse_coordinate(self.coordinates[0], pair = 'left')
            self.parsed_right_coordinate = self._parse_coordinate(self.coordinates[1], pair = 'right')
            self.left_sequence_data = self._UCSC_request(self.parsed_left_coordinate)
            self.right_sequence_data = self._UCSC_request(self.parsed_right_coordinate)

            self.breakpoint_sequence_template = self._build_breakpoint()
            self.breakpoint_primers = self._design_primers(self.breakpoint_sequence_template)
            print(self.breakpoint_primers)
            print(pd.DataFrame.from_dict(self.breakpoint_primers))
            # TODO handle two coordinates for fusions
            pass


    def _parse_coordinate(self, coord, pair = None):
        """Parse provided genomic coordinate for use in UCSC API request. Start and end positions are caluclated
        from themain.py given coordinates, which is assumed to be the center of the desired amplicon.
        Coordinate must be in the format chr1:234,567,890.
        Pass a negative integer for downstream end position"""

        coordinate = coord.lower().split(sep=":")
        # if a single coordinate has been given, start postion is self.template_sequence_length/2.
        if not pair:
            coordinate[1] = int(coordinate[1]) - round(self.template_sequence_length/2)
            coordinate.append(coordinate[1] + self.template_sequence_length)
        elif pair == 'left':
            coordinate[1] = int(coordinate[1]) + self.template_sequence_length * -1
            coordinate.append(coordinate[1] + self.template_sequence_length)
        elif pair == 'right':
            coordinate.append(int(coordinate[1]) + self.template_sequence_length)

        coordinate.append(self.ref_genome)
        coordinate = dict(zip(['chrom', 'start', 'end', 'genome'], coordinate))
        print(coordinate)
        return coordinate

    def _UCSC_request(self):
        """Requests sequence data from UCSC using given coordinate. Returns json."""
        # coords = parse_coordinate(coordinate)
        # coords['genome'] = 'hg38'
        # print(coords)

        # Handle two coordinates on same chromosome
        # break point coordinates will need 2 API requests - make new function for this.
        if self.end_pos:
            end_chrom, end = self.end_pos.split(":")
            start_chrom, start = self.start_pos.split(":")
            # check both coords are on the same chromosome
            if start_chrom != end_chrom:
                raise ChromMismatch("Start and end positions for non-fusions must be on the same chromsome")

            url = f"https://api.genome.ucsc.edu/getData/sequence?genome={self.ref_genome};chrom={start_chrom};start={start};end={end}"


        else:
            # request for single coordinate
            # modify start and end coordinates to be in center of a sequence with length defined by seq_len
            chrom, start = self.start_pos.split(":")
            # given coordinate will be in the center of the returned sequence
            flanking_len = round(self.seq_len/2)
            end = start + flanking_len
            start = start - flanking_len

            url = f"https://api.genome.ucsc.edu/getData/sequence?genome={self.ref_genome};chrom={chrom};start={start};end={end}"

        response = requests.get(url)

        if response.ok:
            print(response.url)
            return response.json()
        else:
            raise BadRequest(f"API request failed with status {response_status_here}")

    def _design_primers(self, template_sequence):
        """design PCR primers using primer3 with default options"""
        # TODO primer design options?
        primers = primer3.design_primers(
            seq_args={'SEQUENCE_ID': 'test', 'SEQUENCE_TEMPLATE': template_sequence},
            global_args={})
        parsed_primers = self._parse_primer3(primers)
        return parsed_primers

    def _build_breakpoint(self):
        '''concantenate breakpoint sequences'''
        breakpoint_sequence = ''.join([self.left_sequence_data['dna'], self.right_sequence_data['dna']])
        return breakpoint_sequence
    def _design_breakpoint_primers(self):
        """placeholder for possible future fusion stuff"""
        #primers = self._design_primers(self.breakpoint_sequence_template)
        #return primers
        pass

    def _parse_primer3(self, primer3_output):
        """parse primer3 output dictionary to be more manageable."""
        parsed = {}

        for key, value in primer3_output.items():
            pair_id = str('PAIR_' + key.split(sep="_")[2])
            if pair_id not in parsed:
                parsed[pair_id] = {}
                parsed[pair_id][key] = value
            else:
                parsed[pair_id][key] = value
        # TODO separate this into new function. Need to think about dict key names to get more readable output df.
        primer_info = {}
        for key, value in parsed.items():
            if key.split(sep="_")[1] in ['0','1','2','3','4','5']:
                primer_info[key] = parsed[key]

        return primer_info
# exceptions
class ChromMismatch(Exception):
    '''raise this if start and end chroms are different in non-breakpoint situation'''
    pass

class BadRequest(Exception):
    '''raise this if API request fails'''
    pass
def prymer_main():
    parser = argparse.ArgumentParser(description='Design PCR primers for given genomic coordinates')
    parser.add_argument('-c',
                        '--coordinates',
                        help='Genomic coordinates in the format chr1:23456',
                        required=True,
                        nargs='*')
    parser.add_argument('-l',
                        '--template_sequence_length',
                        help='Length of template sequence in bp to use for primer design. Default 500bp',
                        type=int,
                        default=500)
    parser.add_argument('-r',
                        '--reference_genome',
                        help='Reference genome to use. Default hg38',
                        default='hg38')
    args = parser.parse_args()


    primer = Primer()
    return primer



# coord1 = "Chr5:12345678"
# coord2 = "Chr5:12345678", "Chr7:12345678"

#test_breakpoint_primers = prymer_main(args.coordinates)
# print(test_primer.sequence_data)
# print(test_primer.sequence_data['dna'])
# print(test_primer.primers)
# print(test_primer.primers['PAIR_0'])
# print(test_primer.primers['PAIR_1']['PRIMER_LEFT_1_SEQUENCE'])

if __name__ == '__main__':
    prymer_main()
