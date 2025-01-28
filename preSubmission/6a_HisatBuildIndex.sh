#!/bin/bash

##########################################################################################
# Notes
##########################################################################################
# since I use array jobs, build index separately from running hisat2

##########################################################################################
# Modules
##########################################################################################
module load HISAT2/3n-20201216-gompi-2022a

##########################################################################################
# Script
##########################################################################################
cat masked/*.fa.masked > maskedgenome.fa
hisat2-build maskedgenome.fa hisat2/maskedgenome