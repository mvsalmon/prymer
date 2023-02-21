def primer3_parser(primer3_results):
   ''' Parse Primer3 designPrimers output, and sort it into a hierachical
   dictionary structure of primer pairs.

   This method return 2 outputs, the list of primer pairs and a dictionary with
   notes (the explanatory output from Primer3).

   Author: Martin CF Thomsen
   '''
   primer_pairs = {}
   notes = {}
   for k in primer3_results:
      if 'PRIMER_RIGHT' == k[:12]:
         key = 'right'
         tmp = k[13:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               primer_pairs[id][key][key2] = primer3_results[k]
            else:
               primer_pairs[id][key]['position'] = primer3_results[k][0]
               primer_pairs[id][key]['length'] = primer3_results[k][1]
         elif tmp[0] == 'EXPLAIN':
            notes[key] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            print(k)
      elif 'PRIMER_LEFT' == k[:11]:
         key = 'left'
         tmp = k[12:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               primer_pairs[id][key][key2] = primer3_results[k]
            else:
               primer_pairs[id][key]['position'] = primer3_results[k][0]
               primer_pairs[id][key]['length'] = primer3_results[k][1]
         elif tmp[0] == 'EXPLAIN':
            notes[key] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            print(k)
      elif 'PRIMER_PAIR' == k[:11]:
         key = 'pair'
         tmp = k[12:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               primer_pairs[id][key][key2] = primer3_results[k]
            else:
               print(k, primer3_results[k])
         elif tmp[0] == 'EXPLAIN':
            notes[key] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            print(k)
      elif 'PRIMER_INTERNAL' == k[:15]:
         key = 'internal'
         tmp = k[16:].split('_', 1)
         if tmp[0].isdigit():
            id = int(tmp[0])
            if not id in primer_pairs:
               primer_pairs[id] = {'pair': {}, 'right': {}, 'left': {},
                                   'internal': {}}
            if len(tmp) > 1:
               key2 = tmp[1].lower()
               primer_pairs[id][key][key2] = primer3_results[k]
            else:
               primer_pairs[id][key]['position'] = primer3_results[k][0]
               primer_pairs[id][key]['length'] = primer3_results[k][1]
         elif tmp[0] == 'EXPLAIN':
            notes['pair'] = primer3_results[k]
         elif tmp == ['NUM','RETURNED']: pass
         else:
            print(k, tmp[0])
      else:
         print(k)

   return list(map(primer_pairs.get, sorted(primer_pairs.keys()))), notes
