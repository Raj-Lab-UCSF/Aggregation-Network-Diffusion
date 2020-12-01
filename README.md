# Aggregation-Network-Diffusion
Code for the Aggregation-Network Diffusion (AND) model of pathology ramification in human brain connectome. There are two main scripts, corresponding to two manuscripts under review.

1) Network Diffusion model jointly for both amyloid and tau, with cross-species interaction terms between the two. This manuscript is currently under review; once it is in print we will update this document with a paper link and citation. The main script that runs this model is 
TauAmyloid_Joint_eNDM_Git.m

We have included group regional amyloid, atrophy and tau SUVr tables for use in this code, however we are unable to post individual subject data. Connectomes and inter-regional distance matrices are also containe dwihtin the folder. All dependencies are resolved, and all files neede to run this scritp are contained within this folder. 

Important: Please add all subfolders within this folder to matlab path prior to running this code.


2) Aggregation and Network Diffusion combined, for tau propagation only

Please cite the paper as follows:
Combined Model of Aggregation And Network Diffusion Recapitulates Alzheimer’s Regional Tau-PET 
Ashish Raj*, Veronica Tora, Xiao Gao, Hanna Cho, Jae Yong Choi, Young Hoon Ryu, Chul Hyoung Lyoo, Bruno Franchi
*Department of Radiology and Biomedical Imaging, University of California at San Francisco

The main script is: franchi_raj_cleaned_Git.m

This version contains code that can do both tau and Abeta
For now we are enabling only dynamics of tau, no Ab (will add later)
This code will reproduce the results of the fully realised Franchi-Raj AND model
We are posting group regional atrophy and tau SUVr tables for use in this code, however we are unable to post individual subject data

Important: Please add all subfolders within this folder to matlab path prior to running this code
