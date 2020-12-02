# Aggregation-Network-Diffusion
Code for the Aggregation-Network Diffusion (AND) model of pathology ramification in human brain connectome. There are two main scripts, corresponding to two manuscripts under review. ALl data, code and other files neededf to rin this code are contained in this folder. There are no other dependencies. The code was tested on MATLAB R2019b on a Windows 10 OS. If using another OS, there should be no problem due to MATLAB's platform-free design. If you want ot run the code on your own data (regional PET or connectomes etc), just replace the variables within the .mat file. The data-bearing variables are named appropriately and intuitively, and are supported by thorough docuemntation within the code. Input/output requirements are listed at teh top of the scripts wihtin the commenting and docuemntation.

Licence: This code is being made freely available under the BSD Clause-3 licence. Please cite related paper when you use the code. 

Installation time should be < 1 min.
Execution: run time shold be  < 1 min (no 3D graphics) or ~3 mins (with 3D graphics enabled - see docuemntation in code). 

1) Network Diffusion model jointly for both amyloid and tau, with cross-species interaction terms between the two. This manuscript is currently under review; once it is in print we will update this document with a paper link and citation. The main script that runs this model is : 
TauAmyloid_Joint_eNDM_Git.m

We have included group regional amyloid, atrophy and tau SUVr tables for use in this code, however we are unable to post individual subject data. Connectomes and inter-regional distance matrices are also contained wihtin the folder. All dependencies are resolved, and all files needed to run this script are contained within this folder. The code should be able to produce any figure conatined in the submited manuscript. 

Important: Please add all subfolders within this folder to matlab path prior to running this code.


2) Aggregation and Network Diffusion combined, for tau propagation only

Please cite the paper as follows:
Combined Model of Aggregation And Network Diffusion Recapitulates Alzheimer’s Regional Tau-PET 
Ashish Raj*, Veronica Tora, Xiao Gao, Hanna Cho, Jae Yong Choi, Young Hoon Ryu, Chul Hyoung Lyoo, Bruno Franchi, biorxiv
*Department of Radiology and Biomedical Imaging, University of California at San Francisco

The main script is: franchi_raj_cleaned_Git.m

This version contains code that can do both tau and Abeta
For now we are enabling only dynamics of tau, no Ab (will add later)
This code will reproduce the results of the fully realised Franchi-Raj AND model
We are posting group regional atrophy and tau SUVr tables for use in this code, however we are unable to post individual subject data

Important: Please add all subfolders within this folder to matlab path prior to running this code
