# TODO command line argument handling
# TODO handle multiple coordinates for fusions
# TODO error handling
from primer import Primer

coord1 = ("Chr5:12345678",)
coord2 = ("Chr5:12345678", "Chr7:12345678")

test_primer = Primer(coord1)
test_breakpoint_primers = Primer(coord2)

print(test_primer.sequence_data)
print(test_primer.sequence_data['dna'])
print(test_primer.primers)
print(test_primer.primers['PAIR_0'])
print(test_primer.primers['PAIR_1']['PRIMER_LEFT_1_SEQUENCE'])


