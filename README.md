# material-composition-lca
Estimating material composition of products using LCI databases
L. Amatuni. 2021
CML, Leiden University

General description:
This software allows estimating the material composition of various products reported in the ecoinvent LCI database.
The software originated withing a scientific work described in the following article where the methodological steps are described in detail: {article submitted to publisher}

How to:
1. Use either Activity Browser or Brightway to create a new project and get access to the ecoinvent database version 3.6
2. Run the appropriate functions from main.py to obtain estimates for a material and product of interest (see the comments for details)
3. If adjustment of the list of plastics is required, explore and run the dict_gen/mat_dict_gen.py file before step 2 

Activity Browser:
A modified version of Activity Browser is provided under the 'add-diffusion' branch of the original repository:
https://github.com/LCA-ActivityBrowser/activity-browser.git
