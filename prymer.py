# TODO command line argument handling
# TODO handle multiple coordinates for fusions
# TODO error handling
import argparse

from primerClass import Primer

def prymer_main():
    parser = argparse.ArgumentParser(description='Design PCR primers for given genomic coordinates')
    parser.add_argument('-c',
                        '--coordinates',
                        help='STRING. Genomic coordinates in the format chr1:23456',
                        required=True,
                        nargs='*')
    parser.add_argument('-l',
                        '--template_sequence_length',
                        help='INT. Length of template sequence in bp to use for primer design. Default = 500',
                        type=int,
                        default=500)
    parser.add_argument('-r',
                        '--reference_genome',
                        help='STRING. Reference genome to use. Default = hg38',
                        default='hg38')
    parser.add_argument('-o',
                        '--output_file',
                        help='STRING. Name of output file. Default = prymer.txt',
                        default='prymer.txt')
    args = parser.parse_args()


    primer = Primer(args)
    return primer



coord1 = "Chr5:12345678"
coord2 = "Chr5:12345678", "Chr7:12345678"

#test_breakpoint_primers = prymer_main(args.coordinates)
# print(test_primer.sequence_data)
# print(test_primer.sequence_data['dna'])
# print(test_primer.primers)
# print(test_primer.primers['PAIR_0'])
# print(test_primer.primers['PAIR_1']['PRIMER_LEFT_1_SEQUENCE'])

if __name__ == '__main__':
    prymer_main()
