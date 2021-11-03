# -*- coding: utf-8 -*-

# Copyright 2021 Netherlands eScience Center and CML, Leiden University
# Licensed under the Apache License, version 2.0. See LICENSE for details.

"""
Created on Wed Jul 17 18:52:09 2019

@author: amatunilt
"""

from brightway2 import *
import numpy as np

projects.set_current("test")
#bw2setup()

db = Database("db")
list_dissip = {"factory"}

#dictionary relates meaningful materials with the corresponding/contributing activities (their keys)
materials_dict = {
        "plastic":
            {"PET": [('db', '243c83e4e16546f58c2f301351a0397d')],
             "PP": [('db', 'fcc4443e7ee64e88b89c9c4d93404d8a')],
            },
        "metals":
            {"copper": [('db', 'e454d24c657340cdb9e20733b8c97a62')]
            }
}

def activity_by_name(name, db):
    candidates = [x for x in db if name in x['name']]
    return candidates[0]

def activity_by_key(key): # tuple(db, key) -> activity (dataset)
    return db.get(key[1])

def list_techno_inventory(lca):
    for k in lca.activity_dict:
        print(activity_by_key(k)["name"], ": ", lca.supply_array[lca.activity_dict[k]])
    print('\n')
    
def composition(materials_dict, lca): #presents product composition from the exisitng lca object based on the metrials dicionary
    for material_group in materials_dict:
        gr_sum = 0
        for material in materials_dict[material_group]:
            mat_sum = 0
            for act_key in materials_dict[material_group][material]:
                mat_sum += lca.supply_array[lca.activity_dict[act_key]]
            print('\t', material, " = ", mat_sum)
            gr_sum += mat_sum
        print(material_group, " = ", gr_sum)
    print('\n')


act = activity_by_name('PC', db)

#assigning 'dissipation' parameter
for exc in act.technosphere():
    exc_name = activity_by_key(exc["input"])["name"]
    exc['dissip'] = 1 if exc_name in list_dissip else 0
    exc.save()
    
#LCI
functional_unit = {act: 1}
method_key = [x for x in methods if "non-renewable resources, copper" in x][0]    
lca = LCA(functional_unit, method_key)

lca.lci()
print("LCI (before): ", list_techno_inventory(lca))
print(composition(materials_dict, lca))


for exc in act.technosphere():
    if  exc['dissip']:
        row = lca.activity_dict[exc["input"]]
        col = lca.activity_dict[act.key]
        lca.technosphere_matrix[row, col] = 0


#output techn. matrix
#for i in range(7):
#    print("-")
#    for j in range(7):
#        print(lca.technosphere_matrix[i, j], end = ' ')

lca.lci_calculation()
lca.lcia()

print("LCI (after): ", list_techno_inventory(lca))
print(method_key, ' : ', lca.score)
print(composition(materials_dict, lca))


