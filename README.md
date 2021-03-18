# covid19-variant-N501Y

## Table of Contents
1. Overview
2. Software Requirements
3. Instructions for Use

## Overview

Forked from mrc-ide/covid19-variant-N501Y:
Originally for paper [Transmission of SARS-CoV-2 Lineage B.1.1.7 in England: Insights from linking epidemiological and genetic data](https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-42-sars-cov-2-variant/)

Now with updated analyses to ask questions about environmental drivers of transmission of Lineage B.1.1.7, for our pre-print [Environmental drivers of SARS-CoV-2 lineage B.1.1.7 transmission in England, October to December 2020](https://doi.org/10.1101/2021.03.09.21253242)

## Software Requirements

Code to run analyses and produce the manuscript figures was written in R version 3.6.3. The analysis code also has some package requirements:

  * tidyverse (1.3.0)
  * lme4 (1.1-23)
  * xtable (1.8-4)

Code has again not been tested with other versions of these packages.

## Instructions for Use

To run all the analyses and generate all figures presented in the paper:

    > Rscript src/temperature_analysis_figures.R
    
To have the latex-formatted results tables output to file, try:

    > Rscript src/temperature_analysis_figures.R > results/results_tables.txt
