# Estimating material composition of products using LCI databases

Developed at the Institute of Environmental Sciences (CML) at Leiden University

## General description:

This software allows estimating the product material composition (PMC) of various goods using ecoinvent LCI database.

The software originated within a scientific work described in the following article where the methodological steps are described in detail:
Amatuni, L., Steubing, B., Heijungs, R., Yamamoto, T., & Mogollón, J. M. (2024). Deriving material composition of products using life cycle inventory databases. Journal of Industrial Ecology, 28(5), 1060-1072. https://doi.org/10.1111/jiec.13538
We highly recommend visiting that work first for a solid understanding and implementation of this product composition estimation algorithm.

## Requirements:

- Running Brightway package in your Python environment
- Access to ecoinvent database (or any other LCI database in Brightway)

## How to:

1. The `pmc algorithm - general.py` file describes the general algorithm for any kind of LCI database you have access to in Brightway environment. Has very detailed comments and explanations. Still, reading the paper above first is highly recommended. Also, the following file might be more applicable if you want to jump straight to real cases and tasks of yours using ecoinvent.
2. The `pmc algorithm - ecoinvent.ipynb` Jupyter Notebook allows you to run all the steps of the algorithm one by one while estimating the content of selected materials in any chosen product in your ecoinvent database in an accessible (for non-experts) format.
3. The `pmc algorithm - ecoinvent.py` file assumes an understanding of the general algorithm and allows you to run the script immediately to estimate specific materials' content in any chosen product in your ecoinvent database. Minimum comments provided.

We are always happy to receive feedback or provide support for the application of our algorithm to the real-world products you work with.
Please, contact: l.t.amatuni@cml.leidenuniv.nl
