Example Dataset:

This is a simulated genotypic-dataset. It includes 1000 genotypes and 200000 variants. The example dataset can be downloaded from DOI: 10.6084/m9.figshare.24025728, Link: https://figshare.com/s/7193d332894f11deb1a5.


A.
For analysis based on additive effects, the user needs to copy MAIN.bed, MAIN.bim, MAIN.fam, and MainPs into the SAFE-h2 folder.

Open the Prog.config for configurations (to select heritability models: the default is EMMAX, see the manual)
run:
python3 SAFE-h2.py


B.
For analysis based on additive and non-additive effects, the user needs to copy MAIN.bed, MAIN.bim, and MAIN.fam (NOT MainPs) into the SAFE-h2 folder.

Open the Prog.config for configurations (to select heritability models: the default is EMMAX, to use user provided association test: the default is plinks glm function, see the manual)
run:
1)   python3 SAFE-h2-preADOH.py                
2)   python3 SAFE-h2-ADOH.pyv

   
C.    
For analysis by minimizing the false-positive contributions, the user needs to copy MAIN.bed, MAIN.bim, MAIN.fam, MainPs, MAIN1.fam, Main1Ps, MAIN2.fam, Main2Ps, MAIN3.fam, Main3Ps, MAIN4.fam, Main4Ps, MAIN5.fam, Main5Ps, MAIN6.fam, and Main6Ps into the SAFE-h2 folder.

Open the Prog.config for configurations (set PURE_Effects (X7) to 1, to select heritability models: the default is EMMAX, see the manual)
run:
python3 SAFE-h2.py





In all of the mentioned analyses, user can also copy the Covar files to include the covariates into the models.





