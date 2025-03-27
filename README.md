# Scripts of analyse and figure generation for the *Phoca vitulina mellonae* environmental DNA study

__Main author:__  Marion Chevrinais  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        Laboratory of genomics   
__Location:__     Maurice Lamontagne Institute, Mont-Joli, Québec, Canada  
__Affiliated publication:__ Chevrinais et al. submitted. Refining the distribution and confirming the phylogenetic relationship of the endangered landlocked harbour seal Phoca vitulina mellonae using environmental DNA. submitted to CJFAS.       
__Contact:__      e-mail: marion.chevrinais@dfo-mpo.gc.ca 

- [Description](#description)
- [Map](#map)
- [Effect of the distance between the sample and seal observations](#effect-of-the-distance-between-the-sample-and-seal-observations)
- [Occupancy modelling](#occupancy-modelling)
- [Phylogenetic tree](#phylogenetic-tree)
- [References](#references)

## Description 

Welcome to this repository with the code necessary to make the analyses and figures of *P. v. mellonae* environmental DNA paper. 

## Map

The code **01_Map.R** to produce maps is in the folder **01_Code**.

Shapefiles to construct the map come from various sources: 
- countries.gpkg = GeoPackage Data
- Waterbody = Ministère des Ressources Naturelles et des Forêts du Québec
- Critical and essential habitats = Fisheries and Oceans Canada
- Watersheds = Ministère des Ressources Naturelles et des Forêts du Québec

Coordinates for Phoca sampling stations are in the **Table S1**. 

The maps produced are used in **Figure 1** and **Figure 2** of the associated manuscript. 

## Effect of the distance between the sample and seal observations

The code **02_Distance_effect.R** to produce **Figure 3**, **Table S9** and model for the effect of distance is in the folder **01_Code**.

## Occupancy modelling

The code **03_Occupancy_modelling.R** to produce **Figure 4**, **Table S10** and model for the occupancy modelling at the sample level is in the folder **01_Code**.
The package *msocc* developped by Stratton et al. 2020 is used for occupancy modelling. 

## Phylogenetic tree

The code **04_Phylo_tree_construction.R** to produce **Figure 5**, **Tables S3 and S4** is in the folder **01_Code**. A maximuum likelihood tree is produced.  

## References

Stratton, C., Sepulveda, A., & Hoegh, A. (2020). msocc: Fit and analyze computationally efficient multi‐scale occupancy models in R. Methods in Ecology and Evolution.
