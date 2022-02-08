# HUVEC

This repo is a collection of scripts that will eventually become an analysis pipeline for the HUVEC data analysis. 

The experiment contains a proliferation phase, wounding phase and thrombin phase. 

Curently, we have written a peak extraction and metrics extraction scripts. 

1) Data is preprocessed using Normalise code, which takes in the raw EXIS outputs, syncronyses the time scale and normalises it to the internal controls using different normalisation methods - min-max normalisation by adding in the controls, Fold Change normalisation, by calculating the average controls starting value and using this as a respective value of FC calculation.

2) peak_extraction.py code takes in the normalised spectral data for the Wounding and Thrombin and extract the peaks and the corresponding metrics. 

