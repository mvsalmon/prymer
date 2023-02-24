# TODO command line argument handling
# TODO handle multiple coordinates for fusions
from primer import Primer

coord1 = ("Chr5:123456",)

test_primer = Primer(coord1)

print(test_primer.sequence_data)
print(test_primer.sequence_data['dna'])
print(test_primer.primers)
print(test_primer.primers['PAIR_0'])
print(test_primer.primers['PAIR_1']['PRIMER'])


