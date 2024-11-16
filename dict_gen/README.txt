The `mat_dict_gen.py` reads the list of all the activities in ecoinevent (their names) in `ecoinvent-3.10-cutoff.csv` 
and outputs their linkage to a specific material group (e.g. plastics) in the `plastics_dict_ecoinvent-3.10-cutoff_m.json` file - the so-called material 
dictionary in our paper.

The 'ecoinvent-3.10-cutoff.csv' file, in turn, is obtained running the db_to_csv(db) function from the 'pmc algorithm - ecoinvent.py' file. Whichever database version you have, the function will output the appropriate list of activities in a CSV format that you can then import into the mat_dict_gen.py for automated processing. 