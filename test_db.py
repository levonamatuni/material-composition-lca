# -*- coding: utf-8 -*-

# Copyright 2021 Netherlands eScience Center and CML, Leiden University
# Licensed under the Apache License, version 2.0. See LICENSE for details.

"""
Created on Wed Jul 17 18:52:09 2019

@author: amatunilt
"""

#from brightway2 import *
import brightway2 as bw
import numpy as np

#CONST
PROD = 'laptop'


bw.projects.set_current("default")
#bw2setup()

db = bw.Database("db")
list_dissip = {"factory"}

#dictionary relates meaningful materials with the corresponding/contributing activities (their keys)
materials_dict = {
        "plastics":
            {"PET": [],
             "PP": [],
            },
        "metals":
            {"copper": [('db', '9d4a7be1e15944dcb936132aea869546')]
            }
}

def activity_by_name(name, db):
    candidates = [x for x in db if name in x['name']]
    return candidates[0]

#obtain activity object given its key (tuple(db, key))
def activity_by_key(key): # tuple(db, key) -> activity (dataset)
    return db.get(key[1])

#list all intermediate (technosphere) flows (activities) in the resulting supply array based on LCI analysis
def list_techno_inventory(lca):
    print("Supply array: ")
    for k in lca.activity_dict:
        print(activity_by_key(k)["name"], ": ", lca.supply_array[lca.activity_dict[k]])
    print('\n')
    
#aggregate and print material fllows from the resulting supply_array based on the predefined 'material dictionary'
def composition(materials_dict, lca): #presents product composition from the exisitng lca object based on the metrials dicionary
    print("Material composition (based on supply array): ")
    for material_group in materials_dict:
        gr_sum = 0
        for material in materials_dict[material_group]:
            mat_sum = 0
            for act_key in materials_dict[material_group][material]:
                mat_sum += lca.supply_array[lca.activity_dict[act_key]]
            print('\t', material, " : ", mat_sum, ' kg')
            gr_sum += mat_sum
        print(material_group, "in total : ", gr_sum, ' kg')
    print('\n')
    
#print technospheere matrix (A)
def print_techno_matrix(lca):
    matrix_size = np.shape(lca.technosphere_matrix)[0]
    for i in range(matrix_size):
        print("\n")
        for j in range(matrix_size):
            print(lca.technosphere_matrix[i, j], end = ' ')         
        
#MAIN code

act = activity_by_name(PROD, db)

#assigning 'dissipation' parameter
for exc in act.technosphere():
    exc_name = activity_by_key(exc["input"])["name"]
    exc['dissip'] = 1 if exc_name in list_dissip else 0
    exc.save()
    
#LCI prior to changing the technosphere matrix (filtering out non-incorporated flows)
print("BEFORE filtering (MC does not make sense):\n")
functional_unit = {act: 1}
method_key = [x for x in bw.methods if "copper" in x][0] #pick the environmental pressure (extracted material/resource from the list of impact methods)
lca = bw.LCA(functional_unit, method_key) 

lca.lci()
lca.lcia()

list_techno_inventory(lca)
print("Material composition (based on inventory vector): \n", method_key, ' : ', lca.score, '\n')
composition(materials_dict, lca)

#LCI afyer changing the technosphere matrix (filtering out non-incorporated flows)
print("AFTER filtering:\n")

#removing the non-incorporated inputs from the technosphere matrix 
for exc in act.technosphere():
    if  exc['dissip']:
        row = lca.activity_dict[exc["input"]]
        col = lca.activity_dict[act.key]
        lca.technosphere_matrix[row, col] = 0

lca.lci_calculation()
lca.lcia()

list_techno_inventory(lca)
print("Material composition (based on inventory vector): \n", method_key, ' : ', lca.score, '\n')
composition(materials_dict, lca)


