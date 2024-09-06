# -*- coding: utf-8 -*-

# Copyright 2024 Netherlands eScience Center and CML, Leiden University
# Licensed under the Apache License, version 2.0. See LICENSE for details.

"""
Created on Thu Aug 22 18:01:26 2019
@author: amatunilt
"""

#NOTES:
#   This script allows to estimate the material composition (MC) of products using real product datasets in ecoinvent database.
#   See paper: Amatuni, L., Steubing, B., Heijungs, R., Yamamoto, T., & Mogollón, J. M. Deriving material composition of products using life cycle inventory databases. Journal of Industrial Ecology. https://doi.org/10.1111/jiec.13538
#   We highly recommend visiting this paper first for a solid understanding of the implementation of this product composition estimation algorithm.
#   We also suggest first reviewing a simpler and more universal script "pmc algorithm - general.py" for detailed comments for each code line.
#   This script is an extension of that general script tailoring it to perform on the ecoinvent database. 
#   The comments under this script will be rather limited. 

#REQUIREMENTS: 
#   To run this script in Spyder or VS Code we first prepared a separate conda environment with an installed Brightway 2 (not 2.5!) framework and the Activity Browser GUI on top.
#   Installing Activity Browser in a separate conda environment as instructed on the AB website will install Brightway and all the packages needed, so this is a good starting point.  
#   Basic understanding of Pythin and Brightway is required to be able to replicate this code for your own products/materials/LCI database  
#   Reading the paper: Amatuni, L., Steubing, B., Heijungs, R., Yamamoto, T., & Mogollón, J. M. Deriving material composition of products using life cycle inventory databases. Journal of Industrial Ecology. https://doi.org/10.1111/jiec.13538
#   Reviewing explanatory comments under the general script first.

#IMPORTS:
from brightway2 import *
from bw2analyzer.matrix_grapher import SparseMatrixGrapher
from functools import cmp_to_key
import matplotlib.pyplot as plt
from bw2data.parameters import ActivityParameter, DatabaseParameter, ProjectParameter, Group
import brightway2 as bw
import numpy as np
import json
import csv
import sys
import os

# Get the directory where the script is located and set the current working directory to the script's directory
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(script_dir)

#CONSTANTS
FLOAT_RND = 5
Cu = "non-renewable resources, copper"
Al = "non-renewable resources, aluminium"
Ta = "non-renewable resources, tantalum"
LAPTOP = 'computer production, laptop'
FRIDGE = 'refrigerator production'
CAR_D = 'passenger car production, diesel'
CAR_P = 'passenger car production, petrol/natural gas'
CAR_E_motor = 'electric motor production, vehicle (electric powertrain)'
SCOOTER_E_motor = 'electric motor production, for electric scooter'
FERRY = 'ferry production'
BUTTER = 'butter production, from cow milk'
db_list = ['cutoff36', 'apos36', 'conseq36']
prod_list = [CAR_P, LAPTOP, FRIDGE]
prod_wght = [1, 3.15, 60] #kg per product from Ecoinvent 3.6
mat_list = [Cu, Al, Ta]
FILE_OUT = 'output.txt'


avoid_activities              = ["treatment", "water", "waste", "container", "box", "packaging", "foam", "electricity", "factory", "adapter", "oxidation", "construction", "heat", "facility", "gas", "freight", "mine", "infrastructure", "conveyor", "road", "building", "used", "maintenance", "transport", "moulding", "mold", "wastewater", "steam", "scrap", "converter"]
avoid_activities_paper          = ["treatment", "water", "waste", "container", "box", "packaging", "foam", "electricity", "factory", "adapter", "oxidation", "construction", "heat", "facility", "gas", "freight", "mine", "infrastructure", "conveyor", "road", "building", "used", "maintenance", "transport", "moulding", "mold", "wastewater", "scrap"]
avoid_activities_scrap        = ["scrap"]
    
materials_dict_cutoff36 = { #this are markets based (since they conveniently include all the regional prod. activities and there's no need to list them separately)
    "metals":
        {"copper_m": [('cutoff36', 'b4f2456cf9cbe7dfeb67c91780bd3e38')],
         "aluminium_cast_alloy_m": [('cutoff36', '9b12c682955e8209a6e7ccdd54afa18e')],
         "aluminium_wrought_alloy_m": [('cutoff36', '88502cc0cdc26c2a9643cd707721b5e4')],
         "aluminium, primary, ingot": [('cutoff36', 'd1297287693145433b8975d85ab1ce82'), ('cutoff36', '97649388d210f81e312b9911734f2d66')],
         "tantalum_m": [('cutoff36', '7aa68196f8d53eb42cdd24dd1725b35c')],
         "cobalt_m": [('cutoff36', '0446f63157ad0192850a65b5ce72bb92')],
         "perm. magnet_m" : [('cutoff36', 'b9233d253b2dbd61f84516076b4ff258')],
         "nyod. ox." : [('cutoff36', '77bfc1edc1346e303862bb8ee401f0fa')]
        }
}

    
materials_dict_apos36 = { #this are markets based (since they conveniently include all the regional prod. activities and there's no need to list them separately)
    "metals":
        {"copper_m": [('apos36', 'b10b9e09d7ed10fa896c067e4c05092b')],
         "aluminium_cast_alloy_m": [('apos36', 'c4e6c1aeba112beeb9540a474e50f78d')],
         "aluminium_wrought_alloy_m": [('apos36', 'c31d45ac03814130c3add51d2c555e63')],
         "aluminium, primary, ingot": [('apos36', 'c7ecee226473eee5b9a7457c9d2bf9b7'), ('apos36', 'ddf457aafeedb02aba959176a13d4269')],
         "tantalum_m": [('apos36', '4ace1aed63e3c81200b0f0ffe8f4de31')]
        }
}
    
materials_dict_conseq36 = { #this are markets based (since they conveniently include all the regional prod. activities and there's no need to list them separately)
    "metals":
        {"copper_m": [('conseq36', '1c4eb36481d69135fefca78435d2c443')],
         "aluminium_cast_alloy_m": [('conseq36', 'f084b0b3e65a60288639cd142803baa6')],
         "aluminium_wrought_alloy_m": [('conseq36', 'dedd8c8513ed33bafe50b54dbfc54140')],
         "aluminium, primary, ingot": [('conseq36', 'dfc6a0f5aa7a87769f2e9969b1aa1415'), ('conseq36', '04367b396c75e63fdb750d3f066d98f6')],
         "tantalum_m": [('conseq36', '0bd909249117a7e7e09f2973cbda6a3a')]
        }
}

materials_dict_toymodel = {
    "plastics": {
        "plastic": [('toy model', '9273337cf4cb4814b60808ddf86a0e55')]
    }
}

#import materials_dict from JSON
with open('dict_gen/plastics_dict_cutoff36_m.json', 'r') as fp:
    materials_dict_cutoff36["plastics"] = json.load(fp) 
    
with open('dict_gen/plastics_dict_apos36_m.json', 'r') as fp:
    materials_dict_apos36["plastics"] = json.load(fp) 
    
with open('dict_gen/plastics_dict_conseq36_m.json', 'r') as fp:
    materials_dict_conseq36["plastics"] = json.load(fp) 
    

def activity_by_name(name, db): #return first activity dataset based on name keyword
    candidates = [x for x in db if name in x['name']]
    candidates = sorted(candidates, key=cmp_to_key(lambda item1, item2: len(item1['name']) - len(item2['name'])))  #shortest name is the best match
    with open(FILE_OUT, 'a') as f:
        print('\n', candidates[0], ' from ', db.name, file=f)
    return candidates[0]

def activity_by_key(key, db): # key = tuple(db, key) -> activity (dataset) in db
    return db.get(key[1])

def list_techno_inventory(lca): #print the resulting total inventory (supply_array) based on LCA
    for k in lca.activity_dict:
        print(activity_by_key(k)["name"], ": ", lca.supply_array[lca.activity_dict[k]])
    print('\n')

def db_inc_filter(db, avoid_activities): #resets inc. parameter in the original db, ~30 min
    i = 0
    prt = -1
    for act in db:
        # update the bar
        i += 1
        pr = int(100*i/len(db))
        if  (pr != prt):
            b = "\rApplying plausible 'unincorporated' parameters to the activities: " + str(pr) + "%"
            print (b, end="\r")
            prt = pr
        for exc in act.technosphere():
            avoid = False
            in_act = bw.get_activity(exc["input"])
            exc_name = in_act._document.product #product name; same as exc["name"] but works for manual db
            # exc_name = activity_by_key(exc["input"], db)["name"] #name of the activity of the exchange
            #if exc["amount"] < 0:
            #    avoid = True
            #else:
            for word in avoid_activities:
                if word in exc_name: #name of the product of the exchange
                    avoid = True
                    break
            exc['incorporated'] = 0.0 if avoid else 1.0
            exc.save()
    print('\n')
    
def db_inc_filter_del_neg(db, avoid_activities): #resets inc. parameter in the original db, ~30 min
    i = 0
    prt = -1
    for act in db:
        # update the bar
        i += 1
        pr = int(100*i/len(db))
        if  (pr != prt):
            b = "\rApplying plausible 'unincorporated' parameters to the activities: " + str(pr) + "%"
            print (b, end="\r")
            prt = pr
        for exc in act.technosphere():
            avoid = False
            in_act = bw.get_activity(exc["input"])
            exc_name = in_act._document.product #product name; same as exc["name"] but works for manual db
            # exc_name = activity_by_key(exc["input"], db)["name"] #name of the activity of the exchange
            if exc["amount"] < 0:
                avoid = True
            else:
                for word in avoid_activities:
                    if word in exc_name: #name of the product of the exchange
                        avoid = True
                        break
            exc['incorporated'] = 0.0 if avoid else 1.0
            exc.save()
    print('\n')
    
def db_inc_reset(db): #all inc -> 1
    i = 0
    prt = -1
    for act in db:
        # update the bar
        i += 1
        pr = int(100*i/len(db))
        if  (pr != prt):
            b = "\rApplying full (1.0) incorporation parameters to the activities: " + str(pr) + "%"
            print (b, end="\r")
            prt = pr
        for exc in act.technosphere():
            exc['incorporated'] = 1.0
            exc.save()
    print('\n')  
            
def lca_exclude_noninc(db, lca): #edit lca techn. matrix (exclude non-incorporative exc) based on incorporation parameter in the db database
    i = 0
    prt = -1
    for act in db:
        # update the bar
        i += 1
        pr = int(100*i/len(db))
        if  (pr != prt):
            b = "\rExcluding the non-incorporated materials from the technosphere matrix: " + str(pr) + "%"
            print (b, end="\r")
            prt = pr
        for exc in act.technosphere():
            try:
                inc = exc['incorporated']
            except:
                inc = 1 #if the incorporation parameter was not entered in ab or set by db_inc_filter() previously
            if  inc < 1:
                row = lca.activity_dict[exc["input"]]
                col = lca.activity_dict[act.key]
                lca.technosphere_matrix[row, col] *= inc
    return lca

def LCA_create(act, FU, material_bioflow, db): #create LCA object based on the reference product in db (production act) and the bioflow of interest     
    functional_unit = {act: FU}
    method_key = [x for x in methods if material_bioflow in x][0]    
    return LCA(functional_unit, method_key)

def composition_bio(act, FU, material_bioflow, db, exclude):
    lca = LCA_create(act, FU, material_bioflow, db)
    lca.lci() #create technosphere
    if exclude:
        lca_exclude_noninc(db, lca) #edit matrix (technosphere)
    lca.lci_calculation() #run LCIA
    lca.lcia()
    with open(FILE_OUT, 'a') as f:
        print('\r', round(lca.score, FLOAT_RND), " kg of ", material_bioflow, " in ", FU, " ", act['name'], end =" ", file=f)
    #SparseMatrixGrapher(lca.biosphere_matrix).graph()
    return lca

def composition_bio_redo(act, FU, db, lca): #redo lcia using same matrices based on act. name in the db (different product)
    lca.redo_lcia({act: FU})
    print('\n', round(lca.score, FLOAT_RND), " kg in ", FU, " ", act['name'])
    
def composition_bio_components(act, FU, material_bioflow, db, exclude):
    lca = composition_bio(act, FU, material_bioflow, db, exclude)
    print('\n', "divided between components:")
    for exc in act.technosphere():
        exc_act = bw.get_activity(exc["input"])
        composition_bio_redo(exc_act, FU * exc["amount"], db, lca)
        #composition_bio(exc_act, FU * exc["amount"], material_bioflow, db, exclude)
    return lca    
    
def composition_tech(act, FU, materials_dict, db, exclude, prod_w): #presents product composition (prod. activity) from the exisitng lca object based on the metrials dicionary
    lca = LCA({act: FU})
    lca.lci() #create technosphere
    if exclude:
        lca_exclude_noninc(db, lca) #edit matrix (technosphere)
    with open(FILE_OUT, 'a') as f:
        print("COMPOSITION\n", file=f)
    lca.lci_calculation() #LCIA calculation
    for material_group in materials_dict:
        print(material_group, ':')
        gr_sum = 0
        for material in materials_dict[material_group]:
            mat_sum = 0
            for act_key in materials_dict[material_group][material]:
                act_key = (act_key[0],act_key[1]) #review!
                mat_sum += lca.supply_array[lca.activity_dict[act_key]]
            with open(FILE_OUT, 'a') as f:
                print(f'  {material:<10}{round(mat_sum, FLOAT_RND):>10} kg OR ', round(mat_sum/prod_w*100, FLOAT_RND), '%', file=f)
            gr_sum += mat_sum
        with open(FILE_OUT, 'a') as f:
            print("\tTOTAL: ", round(gr_sum, FLOAT_RND), ' kg OR ', round(gr_sum/prod_w*100, FLOAT_RND), '%\n', file=f)
    return lca

def db_amount_save(db): #save all the original amounts of the exchanges
    for act in db:
        for exc in act.technosphere():
            exc['amount_save'] = exc['amount']
            exc.save()
            
def db_amount_restore(db): #restore all the original amounts of the exchanges
    for act in db:
        for exc in act.technosphere():
            exc['amount'] = exc['amount_save']
            exc.save()
            
def db_inc_to_amounts(db): #adjust all the amounts of the exchanges based on the incorporation
    for act in db:
        for exc in act.technosphere():
            exc['amount'] = exc['amount_save'] * exc['incorporated']
            exc.save()    
            
# =============================================================================
# def db_inc_reset(db): #restore all the original amounts of the exchanges
#     for act in db:
#         for exc in act.technosphere():
#             try:
#                 t = exc['incorporated']
#             except:
#                 exc['incorporated'] = 1.0
#                 exc.save()
# =============================================================================


def db_to_csv(db): #save datasets (name, key) into the csv file
    #generate a list of act. names and their keys
    list = [['name','key']]
    for act in db:
        list.append([act["name"], act.key[1]])
    #frite into csv
    with open(db.name+'.csv', 'w') as f:
        writer = csv.writer(f, delimiter='|', lineterminator="\n")
        writer.writerows(list)
    return str(db) + " saved into " + db.name + '.csv'

#EXECUTE - START

#bw2setup() #Importing elementary flows, LCIA methods and some other data

projects.set_current("default")
db = Database("consequential310")
db_inc_filter(db, avoid_activities)


act = activity_by_name(LAPTOP, db)

#for exc in act.technosphere():
#   print(bw.get_activity(exc["input"])._document.product, exc['incorporated'])
   
#for exc in act.technosphere():
#   print(exc.as_dict())

# projects.set_current("test1")
# print('\nTest1: Filter out all negative & avoid list\n')
# for prod in prod_list:
#     for db_name in db_list:
#         db = Database(db_name)
#         act = activity_by_name(prod, db)
#         for mat in mat_list:
#             lca = composition_bio(act, 1, mat, db, True)
#             print(' OR ', round(lca.score/prod_wght[prod_list.index(prod)] * 100, FLOAT_RND), '%')

# for prod in prod_list:
#     for db_name in db_list:
#         db = Database(db_name)
#         act = activity_by_name(prod, db)
#         materials_dict = globals()['materials_dict_'+db_name]
#         composition_tech(act, 1, materials_dict, db, True, prod_wght[prod_list.index(prod)])



# projects.set_current('test2') 
# db = Database("cutoff36")
# db_inc_filter(db, avoid_activities_scrap)
# db = Database("apos36")
# db_inc_filter(db, avoid_activities_scrap)
# db = Database("conseq36")
# db_inc_filter(db, avoid_activities_scrap)

# projects.set_current("test2")


# with open(FILE_OUT, 'a') as f:
#     print('\nTest2: Filter out only avoid list\n', file=f)

# for prod in prod_list:
#     for db_name in db_list:
#         db = Database(db_name)
#         act = activity_by_name(prod, db)
#         for mat in mat_list:
#             lca = composition_bio(act, 1, mat, db, True)
#             with open(FILE_OUT, 'a') as f:
#                 print (' OR ', round(lca.score/prod_wght[prod_list.index(prod)] * 100, FLOAT_RND), '%', file=f)
#             #lca = composition_bio(act, 1, mat, db, False)
#             #with open(FILE_OUT, 'a') as f:
#             #    print (' OR ', round(lca.score/prod_wght[prod_list.index(prod)] * 100, FLOAT_RND), '%', file=f)


# for prod in prod_list:
#     for db_name in db_list:
#         db = Database(db_name)
#         act = activity_by_name(prod, db)
#         materials_dict = globals()['materials_dict_'+db_name]
#         composition_tech(act, 1, materials_dict, db, True, prod_wght[prod_list.index(prod)])
#         #composition_tech(act, 1, materials_dict, db, False, prod_wght[prod_list.index(prod)])
        
print('DONE')