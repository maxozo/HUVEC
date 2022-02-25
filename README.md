# HUVEC

![Image](https://github.com/maxozo/HUVEC/blob/main/assets/Experimental_Design.png?raw=true)

This repo is a collection of scripts that will eventually become an analysis pipeline for the HUVEC data analysis. 

The experiment contains a proliferation phase, wounding phase and thrombin phase. 

Curently, we have written a peak extraction and metrics extraction scripts. 

1) Data is preprocessed using Normalise code, which takes in the raw ECIS outputs, syncronyses the time scale and normalises it to the internal controls using different normalisation methods - min-max normalisation by adding in the controls, Fold Change normalisation, by calculating the average controls starting value and using this as a respective value of FC calculation.

2) peak_extraction.py code takes in the normalised spectral data for the Wounding and Thrombin and extract the peaks and the corresponding metrics. 

3) Calculate_Normalisation_Factor.py uses the control measurements to calculate an additional scaling facot for each of the measurements. 

At low frequencies, comes underneath and around the cells. Higher freq - within the cells. Resistence cell layer exposes to the cell with oms law.  
proliferation and barrier function is what we are interested in, ca lover frequency.