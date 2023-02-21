# import requests
# import primer3
#
# class Primer():
#     """Class to manage primer pair objects."""
#
#     def __int__(self):
#         self.left_primer
#         self.right_primer


def parse_primer3(primer3_output):
    """parse primer3 output dictionary"""
    parsed = {}

    for key, value in primer3_output.items():
        pair_id = str('PAIR_' + key.split(sep="_")[2])
        if pair_id not in parsed:
            parsed[pair_id] = {}
            parsed[pair_id][key] = value
        else:
            parsed[pair_id][key] = value

    return parsed