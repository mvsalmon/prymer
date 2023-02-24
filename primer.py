import requests
import primer3


class Primer():
    """Class to manage primer objects.
    seq_len = length of sequence in bp that will be returned from UCSC query. Value will be used to calculate the start
    and end positions of the sequence requested from UCSC.
    When a single coordinate is supplied, this will be at the center of the returned sequence.
    For two coordinates (i.e. fusion), each coordinate should represent the breakpoint, and will be either the
    start or end position of the returned sequence."""

    def __init__(self, coordinates, ref_genome='hg38', seq_len=500):

        if len(coordinates) == 1:
            self.parsed_coordinates = self._parse_coordinate(coordinates[0], ref_genome, seq_len)
        else:
            # TODO handle two coordinates for fusions
            pass

        self.sequence_data = self._UCSC_request()
        self.primers = self._design_primers()

    def _parse_coordinate(self, coord, ref_genome, seq_len):
        """Parse provided genomic coordinate for use in UCSC API request. Start and end positions are caluclated
        from the given coordinates, which is assumed to be the center of the desired amplicon.
        Coordinate must be in the format chr1:234,567,890.
        Pass a negative integer for downstream end position"""

        coordinate = coord.lower().split(sep=":")
        coordinate[1] = int(coordinate[1]) - round(seq_len/2)
        coordinate.append(coordinate[1] + seq_len)
        coordinate.append(ref_genome)
        coordinate = dict(zip(['chrom', 'start', 'end', 'genome'], coordinate))
        print(coordinate)
        return coordinate

    def _UCSC_request(self):
        """Requests sequence data from UCSC using given coordinate. Returns json."""
        # coords = parse_coordinate(coordinate)
        # coords['genome'] = 'hg38'
        # print(coords)
        response = requests.get('https://api.genome.ucsc.edu/getData/sequence', params=self.parsed_coordinates)
        print(response.url)
        return response.json()

    def _design_primers(self):
        """design PCR primers using primer3 with default options"""
        # TODO primer design options?
        primers = primer3.design_primers(
            seq_args={'SEQUENCE_ID': 'test', 'SEQUENCE_TEMPLATE': self.sequence_data['dna']},
            global_args={})
        parsed_primers = self._parse_primer3(primers)
        return parsed_primers

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

        return parsed
