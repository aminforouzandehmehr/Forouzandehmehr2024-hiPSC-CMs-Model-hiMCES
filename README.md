# Forouzandehmehr2024-hiPSC-CMs-Model-hiMCES

The Forouzandehmehr2024 hiPSC-CMs Electro-Mechano-Energetic Model Codes (hiMCES)

**Citation:**

M. Forouzandehmehr, M. Paci, J. Hyttinen, and J. T. Koivumäki, “Mechanisms of hypoxia and contractile dysfunction in ischemia/reperfusion in hiPSC cardiomyocytes: an in silico study,” Disease Models & Mechanisms, p. dmm.050365, Mar. 2024, doi: [10.1242/dmm.050365](https://doi.org/10.1242/dmm.050365).

**Files:**

* `hiMCES.m`: Standard model containing ODEs and state variables
* `hiMCES_SET.m`: SET version used for arrhythmic ischemia-reperfusion simulations
* `MasterCompute_hiMCES.m`: Integrator of the model with modifiable simulation features
* `params_and_states_7_EAD_models.m`: Parameters to generate the SET version of the model
* `rho.m`: rho function (Eq. 15 in the manuscript)
* `BiomarkerCalc.m`: Calls functions to calculate AP, CaT, and Contractile biomarkers (Table S1 in the manuscript)
* `run_hiMCES.m`: This commit will add function to run hiMCES model. It has similar functionality to the MasterCompute_hiMCES, except you can override any single value from the function arguments, without overriding implementation of the function.





