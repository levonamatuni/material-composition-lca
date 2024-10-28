The `mat_dict_gen.py` reads the list of all the activities in ecoinevent (their names) in `ecoinvent-3.10-cutoff.csv` 
and outputs their linkage to a specific material group (e.g. plastics) in the `plastics_dict_ecoinvent-3.10-cutoff_m.json` file - the so-called material 
dictionary in our paper.

- review the scanning function. Could use the one that gives higher score to several keywords matching so that PET is prioritized and assigned over PE for PET keywords 