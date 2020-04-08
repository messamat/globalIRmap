#Intermittent river analysis master script 
#Author: Mathis L. Messager
#Contact info: mathis.messager@mail.mcgill.ca
#Affiliation: 
#Global HydroLAB, Department of Geography, McGill University
#EcoFlows Lab, RiverLy Research Unit, INRAE Lyon

source("R/IRmapping_packages.R") #load packages
source("R/functions.R") #define functions
source("R/plan.R") #creates the drake plan

drake_config(plan, verbose = 2)