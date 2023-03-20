import requests
import primer3
import pandas as pd
import json
class Primer():
    """Class to manage primer objects.
    seq_len = length of sequence in bp that will be returned from UCSC query. Value will be used to calculate the start
    and end positions of the sequence requested from UCSC.
    When a single coordinate is supplied, this will be at the center of the returned sequence.
    For two coordinates (i.e. fusion), each coordinate should represent the breakpoint, and will be either the
    start or end position of the returned sequence."""

    def __init__(self, args):
        self.coordinates = args.coordinates
        self.ref_genome = args.reference_genome
        self.template_sequence_length = args.template_sequence_length
        self.output_file_name = args.output_file

        # process single coordinate
        if len(self.coordinates) == 1:
            self.parsed_coordinate = self._parse_coordinate(self.coordinates[0])
            self.sequence_data = self._UCSC_request(self.parsed_coordinate)
            self.primers = self._design_primers(self.sequence_data['dna'])
            print(self.primers)
            print(pd.DataFrame.from_dict(self.primers))
            self._write_output(self.primers)

        # process pair of coordinates
        elif len(self.coordinates) == 2:
            self.parsed_lef_coordinate = self._parse_coordinate(self.coordinates[0], pair = 'left')
            self.parsed_right_coordinate = self._parse_coordinate(self.coordinates[1], pair = 'right')
            self.left_sequence_data = self._UCSC_request(self.parsed_lef_coordinate)
            self.right_sequence_data = self._UCSC_request(self.parsed_right_coordinate)

            self.breakpoint_sequence_template = self._build_breakpoint()
            self.breakpoint_primers = self._design_breakpoint_primers()
            self._write_output(self.breakpoint_primers)
            # TODO handle two coordinates for fusions
            pass


    def _parse_coordinate(self, coord, pair = None):
        """Parse provided genomic coordinate for use in UCSC API request. Start and end positions are caluclated
        from the given coordinates, which is assumed to be the center of the desired amplicon.
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

    def _UCSC_request(self, coordinate):
        """Requests sequence data from UCSC using given coordinate. Returns json."""
        # coords = parse_coordinate(coordinate)
        # coords['genome'] = 'hg38'
        # print(coords)
        response = requests.get('https://api.genome.ucsc.edu/getData/sequence', params=coordinate)
        print(response.url)
        return response.json()

    def _design_primers(self, template_sequence, sequence_target=[]):
        """design PCR primers using primer3 with default options"""
        # TODO primer design options?
        primers = primer3.design_primers(
            seq_args={'SEQUENCE_ID': 'test',
                      'SEQUENCE_TEMPLATE': template_sequence,
                      'SEQUENCE_TARGET': sequence_target},
            global_args={})
        parsed_primers = self._parse_primer3(primers)
        return parsed_primers

    def _build_breakpoint(self):
        breakpoint_sequence = ''.join([self.left_sequence_data['dna'], self.right_sequence_data['dna']])
        return breakpoint_sequence
    def _design_breakpoint_primers(self):
        """concatenate two sequences either side of a breakpoint to use as template"""
        # define target region in breakpoint sequence to be 10bp either side of break.
        # this could be user defined later.
        target_length = 20
        target_start = self.template_sequence_length - round(target_length/2)
        sequence_target = [target_start, target_length]

        primers = self._design_primers(self.breakpoint_sequence_template, sequence_target)
        return primers

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

    def _write_output(self, primer_data):
        """Save primer3 outpyt to file"""

        with open(self.output_file_name, "w") as out_file:
            json.dump(primer_data, out_file, indent=2)