import requests
import primer3
import primers
coord1 = "Chr1:123456"

def parse_coordinate(coord, end_postition=500):
    """parse provided genomic coordinate for use in UCSC API request. Chromosome and start position are taken from
    the given coordinate. End position (relative to the start position in bp) can also be specified, default = +500bp.
    pass a negative integer for downstream end position"""
    coordinate = coord.lower().split(sep = ":")
    coordinate[1] = int(coordinate[1])
    coordinate.append(coordinate[1] + end_postition)
    coordinate = dict(zip(['chrom', 'start', 'end'], coordinate))
    return coordinate


def UCSC_API_request(coordinate):
    """Request sequence data from UCSC using given coordinate"""
    coords = parse_coordinate(coordinate)
    coords['genome'] = 'hg38'
    print(coords)
    r = requests.get('https://api.genome.ucsc.edu/getData/sequence', params=coords)
    print(r.url)
    return r


query_response = UCSC_API_request(coord1).json()

primers_out = primer3.design_primers(seq_args={'SEQUENCE_ID': 'test', 'SEQUENCE_TEMPLATE': query_response['dna']},
                                 global_args={})

parsed = primers.parse_primer3(primers_out)
