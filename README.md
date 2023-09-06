# Voltage-based-plating-prevention

These Python scripts and related files produce results and figures for the ACS Energy Letters article titled 'Voltage-based strategies for preventing battery degradation under diverse fast-charging conditions' https://doi.org/10.1021/acsenergylett.3c01591. It is intended that the user download the associated datasets from Zenodo https://doi.org/10.5281/zenodo.8298930 and store in the same project folder as the provided files. This code was developed using Spyder v5.4.3 and Python v3.11. A brief description of each script is provided below:

modelCalibration_ExptComparison.py - Figures 1, S1  
sensitivityAnalysisVoltageCorrelation.py - Figures 2, 3, S2, S3  
generateTempCurrentStepParams.py - Algorithm to generate parameter sets with which to run simulations.  
CSVparametersToTxtFile_forCOMSOLParametricSweepImport.py - Convert parameters CSV file into the Txt format required by COMSOL software Parametric Sweep.  
multistepChargeAnalysis_all.py - Figures 4, 5, 6, S4-S7
