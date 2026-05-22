import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.font_manager as font_manager
import os
import os.path
from os import path
import pathlib
import shutil
from pathlib import Path
import pyfiglet
from rich import print
from sys import argv
import time
import pylab



title = pyfiglet.figlet_format('SAFE-h2', width=100)
print(f'[white]{title}[/white]')

print('performs "dynamic Snp Allocation For Estimating narrow-sense Heritability".............')
print('This version is version_2026a.')



Prog1 = pathlib.Path("emmax")
pro11 = "Emmax_model will be used to estimate the heritability."
pro12 = "Emmax_prog not found. To apply Emmax_model, please provide the EMMAX program (binary file) renamed as: emmax"
Prog2 = pathlib.Path("ldak")
pro21 = "Ldak GCTA_models will be used to estimate the heritability."
pro22 = "Ldak Thin_model will be used to estimate the heritability."
pro23 = "Ldak_prog not found. To apply Ldak GCTA and Thin models, please provide the LDAK program (binary file) renamed as: ldak"

Prog3 = pathlib.Path("gcta")
pro31 = "GCTA-GREML_model will be used to estimate the heritability."
pro32 = "GCTA-GREML_prog not found. To apply GCTA-GREML_model, please provide the GCTA-GREML program (binary file) renamed as: gcta"

Prog4 = pathlib.Path("gemma")
pro41 = "GEMMA model will be used to estimate the heritability."
pro42 = "GEMMA_prog not found. To apply GEMMA_model, please provide the GEMMA program (binary file) renamed as: gemma"



Prog1_select = pd.read_csv("Prog.config", sep="\t", usecols=["Impl_status_Emmax", "Impl_status_LdakGCTA", "Impl_status_LdakThin", "Impl_status_GCTA-GREML", "Impl_status_GEMMA", "Test_ADOH", "PURE_Effects"])
Prog1_select_act = Prog1_select['Impl_status_Emmax']
Prog1_select_act.to_csv('Prog1_select_act.txt', header=False, index=False)

Prog2_select = pd.read_csv("Prog.config", sep="\t", usecols=["Impl_status_Emmax", "Impl_status_LdakGCTA", "Impl_status_LdakThin", "Impl_status_GCTA-GREML", "Impl_status_GEMMA", "Test_ADOH", "PURE_Effects"])
Prog2_select_act = Prog2_select['Impl_status_LdakGCTA']
Prog2_select_act.to_csv('Prog2_select_act.txt', header=False, index=False)

Prog3_select = pd.read_csv("Prog.config", sep="\t", usecols=["Impl_status_Emmax", "Impl_status_LdakGCTA", "Impl_status_LdakThin", "Impl_status_GCTA-GREML", "Impl_status_GEMMA", "Test_ADOH", "PURE_Effects"])
Prog3_select_act = Prog3_select['Impl_status_LdakThin']
Prog3_select_act.to_csv('Prog3_select_act.txt', header=False, index=False)

Prog4_select = pd.read_csv("Prog.config", sep="\t", usecols=["Impl_status_Emmax", "Impl_status_LdakGCTA", "Impl_status_LdakThin", "Impl_status_GCTA-GREML", "Impl_status_GEMMA", "Test_ADOH", "PURE_Effects"])
Prog4_select_act = Prog4_select['Impl_status_GCTA-GREML']
Prog4_select_act.to_csv('Prog4_select_act.txt', header=False, index=False)

Prog7_select = pd.read_csv("Prog.config", sep="\t", usecols=["Impl_status_Emmax", "Impl_status_LdakGCTA", "Impl_status_LdakThin", "Impl_status_GCTA-GREML", "Impl_status_GEMMA", "Test_ADOH", "PURE_Effects"])
Prog7_select_act = Prog7_select['Impl_status_GEMMA']
Prog7_select_act.to_csv('Prog7_select_act.txt', header=False, index=False)


Add_OD_Dom = pd.read_csv("Prog.config", sep="\t", usecols=["Impl_status_Emmax", "Impl_status_LdakGCTA", "Impl_status_LdakThin", "Impl_status_GCTA-GREML", "Impl_status_GEMMA", "Test_ADOH", "PURE_Effects"])
Add_OD_Dom_act = Add_OD_Dom['Test_ADOH']
Add_OD_Dom_act.to_csv('Add_OD_Dom.txt', header=False, index=False)


valueA = '1'  
valueAa = '0'
valueB = '1'
valueBb = '0'
valueC = '1'
valueCc = '0'
valueD = '1'  
valueDd = '0'
        
valueE = 'plink_based_all'  
valueEe = 'initial_files'
valueEee = 'intermediate_files'


valueG = '1'
valueGg = '0'





print('SAFE-h2 continues...checking for inputs')

bed_check = pathlib.Path("MAIN.bed")
bed1 = "Main.bed is harvested"
bed2 = "SAFE-h2 stops working: Main.bed not found!"
fam_check = pathlib.Path("MAIN.fam")
fam1 = "Main.fam is harvested"
fam2 = "SAFE-h2 stops working: Main.fam not found!"
bim_check = pathlib.Path("MAIN.bim")
bim1 = "Main.bim is harvested"
bim2 = "SAFE-h2 stops working: Main.bim not found!"


if bed_check.exists ():
    print(bed1)
else:
    print(bed2)
    exit()

if fam_check.exists ():
    print(fam1)
else:
    print(fam2)
    exit()

if bim_check.exists ():
    print(bim1)
else:
    print(bim2)
    exit()




valueE = 'plink_based_all'  
valueEe = 'initial_files'
valueEee = 'intermediate_files'
valueEeee = 'binary_observations'



with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        print(' "SAFE-h2 builds VCFs and bfiles for 2 Dominance, 2 Overdominance, and 1 heterosis scenarios"...')
        
with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[Add] --out Additive; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[Dom] --out Dominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[OD] --out Overdominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[dom] --out dominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[od] --out overdominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[het] --out Heter"
        
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+0/1+1/1+g' Dominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+0/1+0/0+g' dominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+1/1+x/x+g' Overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+0/1+1/1+g' Overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+x/x+0/1+g' Overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+0/0+x/x+g' overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+0/1+0/0+g' overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+x/x+0/1+g' overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "sed -i 's+1/1+0/0+g' Heter.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())



        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "./plink2 --vcf Additive.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Additive; ./plink2 --vcf Dominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Dominance; ./plink2 --vcf Overdominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Overdominance; ./plink2 --vcf dominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out dominance; ./plink2 --vcf overdominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out overdominance; ./plink2 --vcf Heter.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Heter"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "cp -f MAIN.fam Additive.fam; cp -f MAIN.fam Dominance.fam; cp -f MAIN.fam Overdominance.fam; cp -f MAIN.fam dominance.fam; cp -f MAIN.fam overdominance.fam; cp -f MAIN.fam Heter.fam"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())




with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        print('"SAFE-h2 has finalized building VCFs and bfiles for 2 Dominance, 2 Overdominance and 1 heterosis scenarios. The SAFE-h2 version is 2026a."')

        print('SAFE-h2 Citation: Behrooz Darbani, Mogens Nicolaisen. On the genetic origions of phenotypes in genome-wide association studies: the SAFE-h2 tool for exploring additive and non-additive allelic effects. BMC Bioinformatics 2026. https://doi.org/10.1186/s12859-026-06464-6')

        print('PLINK Citation: GigaScience 2015Dec;4(1):s13742-015-0047-8 (doi:https://doi.org/10.1186/s13742-015-0047-8)')

        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "rm Prog1_select_act.txt; rm Prog2_select_act.txt; rm Prog3_select_act.txt; rm Prog4_select_act.txt; rm Dominance.vcf; rm Overdominance.vcf; rm dominance.vcf; rm overdominance.vcf; rm Heter.vcf; rm *.nosex; rm *.tped; rm *.tfam; rm *.map; rm *.BN.kinf; rm *.log; rm Dominance.reml; m Overdominance.reml; rm Additive.reml; rm AD_MergePs1; rm AD_MergePs2; rm AD_MergePs3; rm ADO_MergePs1; rm ADO_MergePs2; rm ADO_MergePs3; rm Additive.vcf; rm AdditiveEx.vcf; rm AdditiveEx1.vcf; rm DominanceEx.vcf; rm DominanceEx1.vcf; rm OverdominanceEx.vcf; rm OverdominanceEx1.vcf; rm Merged_Add-Dom1.vcf; rm Merged_Add-Dom2.vcf; rm Merged_Add-Dom3.vcf; rm Merged_Add-Dom_OD1.vcf; rm Merged_Add-Dom_OD2.vcf; rm Merged_Add-Dom_OD3.vcf; rm headers"
        
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEe in lineE:
        command = "mkdir SAFE-h2_output_combinatory-bfiles; mv Additive* SAFE-h2_output_combinatory-bfiles; mv Dominance* SAFE-h2_output_combinatory-bfiles; mv dominance* SAFE-h2_output_combinatory-bfiles; mv Overdominance* SAFE-h2_output_combinatory-bfiles; mv overdominance* SAFE-h2_output_combinatory-bfiles; mv Heter* SAFE-h2_output_combinatory-bfiles; rm Add_OD_Dom.txt"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        exit()












        
with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        print(' "SAFE-h2 builds VCFs and bfiles for 2 Dominance, 2 Overdominance, and 1 heterosis scenarios"...')
        
with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[Add] --out Additive; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[Dom] --out Dominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[OD] --out Overdominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[dom] --out dominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[od] --out overdominance; ./plink2 --bfile MAIN --recode vcf --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --new-id-max-allele-len 20 truncate -set-all-var-ids @:#\$r:\$a[het] --out Heter"
        
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+0/1+1/1+g' Dominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        
with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+0/1+0/0+g' dominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


        
with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+1/1+x/x+g' Overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+0/1+1/1+g' Overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+x/x+0/1+g' Overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+0/0+x/x+g' overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+0/1+0/0+g' overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+x/x+0/1+g' overdominance.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed -i 's+1/1+0/0+g' Heter.vcf"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "./plink2 --vcf Additive.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Additive; ./plink2 --vcf Dominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Dominance; ./plink2 --vcf Overdominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Overdominance; ./plink2 --vcf dominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out dominance; ./plink2 --vcf overdominance.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out overdominance; ./plink2 --vcf Heter.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out Heter"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "cp -f MAIN.fam Additive.fam; cp -f MAIN.fam Dominance.fam; cp -f MAIN.fam Overdominance.fam; cp -f MAIN.fam dominance.fam; cp -f MAIN.fam overdominance.fam; cp -f MAIN.fam Heter.fam"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        







with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        print(' "SAFE-h2 utilizes Plink2 (--glm) to performe association analysis for Additive, 2 Dominance, 2 Overdominance, and 1 heterosis scenarios and using scenario-specific genotype adjustments"...')



fileCovar_plink = pathlib.Path("Covar_plink")

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        if fileCovar_plink.exists ():
            print('Association analysis by plink (--glm)...covariates included...')
            command = "./plink2 --bfile Additive --allow-no-sex --chr-set 90 --allow-extra-chr --covar Covar_plink --glm --out Additive --seed 1551416198; ./plink2 --bfile Dominance --allow-no-sex --chr-set 90 --allow-extra-chr --covar Covar_plink --glm --out Dominance --seed 1551416198; ./plink2 --bfile Overdominance --allow-no-sex --chr-set 90 --allow-extra-chr --covar Covar_plink --glm --out Overdominance --seed 1551416198; ./plink2 --bfile dominance --allow-no-sex --chr-set 90 --allow-extra-chr --covar Covar_plink --glm --out dominance --seed 1551416198; ./plink2 --bfile overdominance --allow-no-sex --chr-set 90 --allow-extra-chr --covar Covar_plink --glm --out overdominance --seed 1551416198; ./plink2 --bfile Heter --allow-no-sex --chr-set 90 --allow-extra-chr --covar Covar_plink --glm --out Heter --seed 1551416198; gawk -i inplace '!/Cov00/' Additive.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' Dominance.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' Overdominance.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' dominance.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' overdominance.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' Heter.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' Additive.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' Dominance.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' Overdominance.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' dominance.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' overdominance.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' Heter.PHENO1.glm.logistic"
            ret = subprocess.run(command, capture_output=True, shell=True)
            # before Python 3.7:
            # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
            print(ret.stdout.decode())
        else:
            print('Association analysis by plink (--glm)...without covariates...')
            command = "./plink2 --bfile Additive --allow-no-sex --chr-set 90 --allow-extra-chr --glm --out Additive --seed 1551416198; ./plink2 --bfile Dominance --allow-no-sex --chr-set 90 --allow-extra-chr --glm --out Dominance --seed 1551416198; ./plink2 --bfile Overdominance --allow-no-sex --chr-set 90 --allow-extra-chr --glm --out Overdominance --seed 1551416198; ./plink2 --bfile dominance --allow-no-sex --chr-set 90 --allow-extra-chr --glm --out dominance --seed 1551416198; ./plink2 --bfile overdominance --allow-no-sex --chr-set 90 --allow-extra-chr --glm --out overdominance --seed 1551416198; ./plink2 --bfile Heter --allow-no-sex --chr-set 90 --allow-extra-chr --glm --out Heter --seed 1551416198"
            ret = subprocess.run(command, capture_output=True, shell=True)
            # before Python 3.7:
            # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
            print(ret.stdout.decode())



with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "mv Additive.PHENO1.glm.linear Additive.Ps; mv Dominance.PHENO1.glm.linear Dominance.Ps; mv Overdominance.PHENO1.glm.linear Overdominance.Ps; mv dominance.PHENO1.glm.linear dominance.Ps; mv overdominance.PHENO1.glm.linear overdominance.Ps; mv Heter.PHENO1.glm.linear Heter.Ps; mv Additive.PHENO1.glm.logistic Additive.Ps; mv Dominance.PHENO1.glm.logistic Dominance.Ps; mv Overdominance.PHENO1.glm.logistic Overdominance.Ps; mv dominance.PHENO1.glm.logistic dominance.Ps; mv overdominance.PHENO1.glm.logistic overdominance.Ps; mv Heter.PHENO1.glm.logistic Heter.Ps"
        
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "awk '{print $3,$9,$12}' Additive.Ps > AdditivePs0; awk '{print $3,$9,$12}' Dominance.Ps > DominancePs0; awk '{print $3,$9,$12}' Overdominance.Ps > OverdominancePs0; awk '{print $3,$9,$12}' dominance.Ps > dominancePs0; awk '{print $3,$9,$12}' overdominance.Ps > overdominancePs0; awk '{print $3,$9,$12}' Heter.Ps > HeterPs0"
        
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "sed '1d' AdditivePs0 > AdditivePs; sed '1d' DominancePs0 > DominancePs; sed '1d' OverdominancePs0 > OverdominancePs; sed '1d' dominancePs0 > dominancePs; sed '1d' overdominancePs0 > overdominancePs; sed '1d' HeterPs0 > HeterPs"
        
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())





with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        print(' "SAFE-h2 builds intermediate files for Additive+Dominance, Additive+Dominance+Overdominance, and Additive+Dominance+Overdominance+Heterosis scenarios"...')

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "cat AdditivePs DominancePs dominancePs > AD_MergePs; sed -i 's/NA/0.99999999999999/' AD_MergePs; sed -i 's/NA/0.99999999999999/' AD_MergePs"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = """awk 'sub(/\[/," "){$1=$1}1' OFS='\t' AD_MergePs > AD_MergePs1x"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub("e","\t",$4)}1' AD_MergePs1x > AD_MergePs1xx"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub(" ","\t")}1' AD_MergePs1xx > AD_MergePs1"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        AD_MergePs1_check = pathlib.Path("AD_MergePs1")
        columns = ['ID1', 'ID2', 'beta', 'pvalue1', 'pvalue2']
        df0_AD_MergePs1 = pd.read_csv("AD_MergePs1", header=None, sep='\t', names=columns)
        df0_AD_MergePs2 = df0_AD_MergePs1.sort_values(by = ['pvalue2', 'pvalue1', 'ID2'])
        df0_AD_MergePs2.to_csv('AD_MergePs2', header=False, sep='\t', index=False) 



with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = """awk '!seen[$1]++' AD_MergePs2 > AD_MergePs3; awk '{print $1"["$2}' < AD_MergePs3 > AD_MergePs4"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())



with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "cat AdditivePs DominancePs dominancePs OverdominancePs overdominancePs > ADO_MergePs; sed -i 's/NA/0.99999999999999/' ADO_MergePs; sed -i 's/NA/0.99999999999999/' ADO_MergePs"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = """awk 'sub(/\[/," "){$1=$1}1' OFS='\t' ADO_MergePs > ADO_MergePs1x"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub("e","\t",$4)}1' ADO_MergePs1x > ADO_MergePs1xx"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub(" ","\t")}1' ADO_MergePs1xx > ADO_MergePs1"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        ADO_MergePs1_check = pathlib.Path("ADO_MergePs1")
        columns = ['ID1', 'ID2', 'beta', 'pvalue1', 'pvalue2']
        df0_ADO_MergePs1 = pd.read_csv("ADO_MergePs1", header=None, sep='\t', names=columns)
        df0_ADO_MergePs2 = df0_ADO_MergePs1.sort_values(by = ['pvalue2', 'pvalue1', 'ID2'])
        df0_ADO_MergePs2.to_csv('ADO_MergePs2', header=False, sep='\t', index=False) 

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = """awk '!seen[$1]++' ADO_MergePs2 > ADO_MergePs3; awk '{print $1"["$2}' < ADO_MergePs3 > ADO_MergePs4"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())





with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = "cat AdditivePs DominancePs dominancePs OverdominancePs overdominancePs HeterPs > ADOH_MergePs; sed -i 's/NA/0.99999999999999/' ADOH_MergePs; sed -i 's/NA/0.99999999999999/' ADOH_MergePs"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = """awk 'sub(/\[/," "){$1=$1}1' OFS='\t' ADOH_MergePs > ADOH_MergePs1x"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub("e","\t",$4)}1' ADOH_MergePs1x > ADOH_MergePs1xx"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub(" ","\t")}1' ADOH_MergePs1xx > ADOH_MergePs1"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        
with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        ADOH_MergePs1_check = pathlib.Path("ADOH_MergePs1")
        columns = ['ID1', 'ID2', 'beta', 'pvalue1', 'pvalue2']
        df0_ADOH_MergePs1 = pd.read_csv("ADOH_MergePs1", header=None, sep='\t', names=columns)
        df0_ADOH_MergePs2 = df0_ADOH_MergePs1.sort_values(by = ['pvalue2', 'pvalue1', 'ID2'])
        df0_ADOH_MergePs2.to_csv('ADOH_MergePs2', header=False, sep='\t', index=False) 


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        command = """awk '!seen[$1]++' ADOH_MergePs2 > ADOH_MergePs3; awk '{print $1"["$2}' < ADOH_MergePs3 > ADOH_MergePs4"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


















        





with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        print(' "SAFE-h2 builds intermediate files for Additive+Dominance, Additive+Dominance+Overdominance, Additive+Dominance+Overdominance+Heterosis scenarios"...')


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = "mv SAFE-h2_output_combinatory-bfiles/* .; rmdir SAFE-h2_output_combinatory-bfiles"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


        
with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = "cat AdditivePs DominancePs dominancePs > AD_MergePs; sed -i 's/NA/0.99999999999999/' AD_MergePs; sed -i 's/NA/0.99999999999999/' AD_MergePs"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = """awk 'sub(/\[/," "){$1=$1}1' OFS="\t" AD_MergePs > AD_MergePs1x"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub("e","\t",$4)}1' AD_MergePs1x > AD_MergePs1xx"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub(" ","\t")}1' AD_MergePs1xx > AD_MergePs1"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        AD_MergePs1_check = pathlib.Path("AD_MergePs1")
        columns = ['ID1', 'ID2', 'beta', 'pvalue1', 'pvalue2']
        df0_AD_MergePs1 = pd.read_csv("AD_MergePs1", header=None, sep='\t', names=columns)
        df0_AD_MergePs2 = df0_AD_MergePs1.sort_values(by = ['pvalue2', 'pvalue1', 'ID2'])
        df0_AD_MergePs2.to_csv('AD_MergePs2', header=False, sep='\t', index=False) 



with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = """awk '!seen[$1]++' AD_MergePs2 > AD_MergePs3; awk '{print $1"["$2}' < AD_MergePs3 > AD_MergePs4"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())




with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = "cat AdditivePs DominancePs dominancePs OverdominancePs overdominancePs > ADO_MergePs; sed -i 's/NA/0.99999999999999/' ADO_MergePs; sed -i 's/NA/0.99999999999999/' ADO_MergePs"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = """awk 'sub(/\[/," "){$1=$1}1' OFS="\t" ADO_MergePs > ADO_MergePs1x"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub("e","\t",$4)}1' ADO_MergePs1x > ADO_MergePs1xx"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub(" ","\t")}1' ADO_MergePs1xx > ADO_MergePs1"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        ADO_MergePs1_check = pathlib.Path("ADO_MergePs1")
        columns = ['ID1', 'ID2', 'beta', 'pvalue1', 'pvalue2']
        df0_ADO_MergePs1 = pd.read_csv("ADO_MergePs1", header=None, sep='\t', names=columns)
        df0_ADO_MergePs2 = df0_ADO_MergePs1.sort_values(by = ['pvalue2', 'pvalue1', 'ID2'])
        df0_ADO_MergePs2.to_csv('ADO_MergePs2', header=False, sep='\t', index=False) 



with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = """awk '!seen[$1]++' ADO_MergePs2 > ADO_MergePs3; awk '{print $1"["$2}' < ADO_MergePs3 > ADO_MergePs4"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())






with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = "cat AdditivePs DominancePs dominancePs OverdominancePs overdominancePs HeterPs > ADOH_MergePs; sed -i 's/NA/0.99999999999999/' ADOH_MergePs; sed -i 's/NA/0.99999999999999/' ADOH_MergePs"
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = """awk 'sub(/\[/," "){$1=$1}1' OFS="\t" ADOH_MergePs > ADOH_MergePs1x"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub("e","\t",$4)}1' ADOH_MergePs1x > ADOH_MergePs1xx"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())
        command = """awk '{gsub(" ","\t")}1' ADOH_MergePs1xx > ADOH_MergePs1"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())


with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        ADOH_MergePs1_check = pathlib.Path("ADOH_MergePs1")
        columns = ['ID1', 'ID2', 'beta', 'pvalue1', 'pvalue2']
        df0_ADOH_MergePs1 = pd.read_csv("ADOH_MergePs1", header=None, sep='\t', names=columns)
        df0_ADOH_MergePs2 = df0_ADOH_MergePs1.sort_values(by = ['pvalue2', 'pvalue1', 'ID2'])
        df0_ADOH_MergePs2.to_csv('ADOH_MergePs2', header=False, sep='\t', index=False)

with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueEee in lineE:
        command = """awk '!seen[$1]++' ADOH_MergePs2 > ADOH_MergePs3; awk '{print $1"["$2}' < ADOH_MergePs3 > ADOH_MergePs4"""
        ret = subprocess.run(command, capture_output=True, shell=True)
        # before Python 3.7:
        # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        print(ret.stdout.decode())










print(' "SAFE-h2 builds combinatory VCF"...')





command = "./plink2 --bfile Additive --recode vcf --extract AD_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out AdditiveEx_AD; ./plink2 --bfile Dominance --recode vcf --extract AD_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out DominanceEx_AD; ./plink2 --bfile dominance --recode vcf --extract AD_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out dominanceEx_AD"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "gawk -i inplace '!/#/' AdditiveEx_AD.vcf; gawk -i inplace '!/#/' DominanceEx_AD.vcf; gawk -i inplace '!/#/' dominanceEx_AD.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "mv AdditiveEx_AD.vcf AdditiveEx_AD1.vcf; mv DominanceEx_AD.vcf DominanceEx_AD1.vcf; mv dominanceEx_AD.vcf dominanceEx_AD1.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "cat AdditiveEx_AD1.vcf DominanceEx_AD1.vcf dominanceEx_AD1.vcf > Merged_Add-Dom1.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "sort -k1,1V -k2,2n Merged_Add-Dom1.vcf > Merged_Add-Dom3.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "./plink2 --bfile Additive --recode vcf --allow-no-sex --chr-set 90 --allow-extra-chr --out head; grep -E '#' head.vcf > headers; cat headers Merged_Add-Dom3.vcf > Merged_Add-Dom_final.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "./plink2 --vcf Merged_Add-Dom_final.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out MAIN_AD; cp -f Additive.fam MAIN_AD.fam"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())




  

command = "./plink2 --bfile Additive --recode vcf --extract ADO_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out AdditiveEx_ADO; ./plink2 --bfile Dominance --recode vcf --extract ADO_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out DominanceEx_ADO; ./plink2 --bfile Overdominance --recode vcf --extract ADO_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out OverdominanceEx_ADO; ./plink2 --bfile dominance --recode vcf --extract ADO_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out dominanceEx_ADO; ./plink2 --bfile overdominance --recode vcf --extract ADO_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out overdominanceEx_ADO"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "gawk -i inplace '!/#/' AdditiveEx_ADO.vcf; gawk -i inplace '!/#/' DominanceEx_ADO.vcf; gawk -i inplace '!/#/' OverdominanceEx_ADO.vcf; gawk -i inplace '!/#/' dominanceEx_ADO.vcf; gawk -i inplace '!/#/' overdominanceEx_ADO.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "mv AdditiveEx_ADO.vcf AdditiveEx_ADO1.vcf; mv DominanceEx_ADO.vcf DominanceEx_ADO1.vcf; mv OverdominanceEx_ADO.vcf OverdominanceEx_ADO1.vcf; mv dominanceEx_ADO.vcf dominanceEx_ADO1.vcf; mv overdominanceEx_ADO.vcf overdominanceEx_ADO1.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "cat AdditiveEx_ADO1.vcf DominanceEx_ADO1.vcf dominanceEx_ADO1.vcf OverdominanceEx_ADO1.vcf overdominanceEx_ADO1.vcf > Merged_Add-Dom-OD1.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "sort -k1,1V -k2,2n Merged_Add-Dom-OD1.vcf > Merged_Add-Dom-OD3.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "cat headers Merged_Add-Dom-OD3.vcf > Merged_Add-Dom-OD_final.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "./plink2 --vcf Merged_Add-Dom-OD_final.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out MAIN_ADO; cp -f Additive.fam MAIN_ADO.fam"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())










command = "./plink2 --bfile Additive --recode vcf --extract ADOH_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out AdditiveEx_ADOH; ./plink2 --bfile Dominance --recode vcf --extract ADOH_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out DominanceEx_ADOH; ./plink2 --bfile Overdominance --recode vcf --extract ADOH_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out OverdominanceEx_ADOH; ./plink2 --bfile dominance --recode vcf --extract ADOH_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out dominanceEx_ADOH; ./plink2 --bfile overdominance --recode vcf --extract ADOH_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out overdominanceEx_ADOH; ./plink2 --bfile Heter --recode vcf --extract ADOH_MergePs4 --allow-no-sex --chr-set 90 --allow-extra-chr --out HeterEx_ADOH"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "gawk -i inplace '!/#/' AdditiveEx_ADOH.vcf; gawk -i inplace '!/#/' DominanceEx_ADOH.vcf; gawk -i inplace '!/#/' OverdominanceEx_ADOH.vcf; gawk -i inplace '!/#/' dominanceEx_ADOH.vcf; gawk -i inplace '!/#/' overdominanceEx_ADOH.vcf; gawk -i inplace '!/#/' HeterEx_ADOH.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "mv AdditiveEx_ADOH.vcf AdditiveEx_ADOH1.vcf; mv DominanceEx_ADOH.vcf DominanceEx_ADOH1.vcf; mv OverdominanceEx_ADOH.vcf OverdominanceEx_ADOH1.vcf; mv dominanceEx_ADOH.vcf dominanceEx_ADOH1.vcf; mv overdominanceEx_ADOH.vcf overdominanceEx_ADOH1.vcf; mv HeterEx_ADOH.vcf HeterEx_ADOH1.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "cat AdditiveEx_ADOH1.vcf DominanceEx_ADOH1.vcf dominanceEx_ADOH1.vcf OverdominanceEx_ADOH1.vcf overdominanceEx_ADOH1.vcf HeterEx_ADOH1.vcf > Merged_Add-Dom-OD-Het1.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "sort -k1,1V -k2,2n Merged_Add-Dom-OD-Het1.vcf > Merged_Add-Dom-OD-Het3.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "cat headers Merged_Add-Dom-OD-Het3.vcf > Merged_Add-Dom-OD-Het_final.vcf"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())

command = "./plink2 --vcf Merged_Add-Dom-OD-Het_final.vcf --make-bed --allow-no-sex --chr-set 90 --allow-extra-chr --out MAIN_ADOH; cp -f Additive.fam MAIN_ADOH.fam"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())






with open("Add_OD_Dom.txt", "a+") as ActE:
    ActE.seek(0) # set position to start of file
    lineE = ActE.read().splitlines()
    if valueE in lineE:
        if fileCovar_plink.exists ():
            print(' "Association analysis by plink (--glm)...covariates included...on the final combinatory genotypic data"...')
            command = "./plink2 --bfile MAIN --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --glm --covar Covar_plink --seed 1551416198 --out MAIN; ./plink2 --bfile MAIN_AD --allow-no-sex --chr-set 90 --allow-extra-chr --glm --covar Covar_plink --seed 1551416198 --out MAIN_AD; ./plink2 --bfile MAIN_ADO --allow-no-sex --chr-set 90 --allow-extra-chr --glm --covar Covar_plink --seed 1551416198 --out MAIN_ADO; ./plink2 --bfile MAIN_ADOH --allow-no-sex --chr-set 90 --allow-extra-chr --glm --covar Covar_plink --seed 1551416198 --out MAIN_ADOH; gawk -i inplace '!/Cov00/' MAIN.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' MAIN_AD.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' MAIN_ADO.PHENO1.glm.linear; gawk -i inplace '!/Cov00/' MAIN_ADOH.PHENO1.glm.linear; mv MAIN.PHENO1.glm.linear MAIN.Ps; mv MAIN_AD.PHENO1.glm.linear MAIN_AD.Ps; mv MAIN_ADO.PHENO1.glm.linear MAIN_ADO.Ps; mv MAIN_ADOH.PHENO1.glm.linear MAIN_ADOH.Ps; gawk -i inplace '!/Cov00/' MAIN.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' MAIN_AD.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' MAIN_ADO.PHENO1.glm.logistic; gawk -i inplace '!/Cov00/' MAIN_ADOH.PHENO1.glm.logistic; mv MAIN.PHENO1.glm.logistic MAIN.Ps; mv MAIN_AD.PHENO1.glm.logistic MAIN_AD.Ps; mv MAIN_ADO.PHENO1.glm.logistic MAIN_ADO.Ps; mv MAIN_ADOH.PHENO1.glm.logistic MAIN_ADOH.Ps; awk '{print $3,$9,$12}' MAIN.Ps > MAINPs0; awk '{print $3,$9,$12}' MAIN_AD.Ps > MAIN_ADPs0; awk '{print $3,$9,$12}' MAIN_ADO.Ps > MAIN_ADOPs0; awk '{print $3,$9,$12}' MAIN_ADOH.Ps > MAIN_ADOHPs0; sed -i 's/NA/0.99999999999999/' MAINPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOHPs0; sed -i 's/NA/0.99999999999999/' MAINPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOHPs0; sed '1d' MAINPs0 > MainPs_Add1; sed '1d' MAIN_ADPs0 > MainPs_AddDom1; sed '1d' MAIN_ADOPs0 > MainPs_AddDomOD1; sed '1d' MAIN_ADOHPs0 > MainPs_AddDomODHet1; sed -e 's/ /\t/g' MainPs_Add1 > MainPs_Add; sed -e 's/ /\t/g' MainPs_AddDom1 > MainPs_AddDom; sed -e 's/ /\t/g' MainPs_AddDomOD1 > MainPs_AddDomOD; sed -e 's/ /\t/g' MainPs_AddDomODHet1 > MainPs_AddDomODHet"
            
            ret = subprocess.run(command, capture_output=True, shell=True)
            # before Python 3.7:
            # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
            print(ret.stdout.decode())
        else:
            print(' "Association analysis by plink (--glm)...without covariates...on the final combinatory genotypic data"...')
            command = "./plink2 --bfile MAIN --allow-no-sex --max-alleles 2 --chr-set 90 --allow-extra-chr --glm --seed 1551416198 --out MAIN; ./plink2 --bfile MAIN_AD --allow-no-sex --chr-set 90 --allow-extra-chr --glm --seed 1551416198 --out MAIN_AD; ./plink2 --bfile MAIN_ADO --allow-no-sex --chr-set 90 --allow-extra-chr --glm --seed 1551416198 --out MAIN_ADO; ./plink2 --bfile MAIN_ADOH --allow-no-sex --chr-set 90 --allow-extra-chr --glm --seed 1551416198 --out MAIN_ADOH; mv MAIN.PHENO1.glm.linear MAIN.Ps; mv MAIN_AD.PHENO1.glm.linear MAIN_AD.Ps; mv MAIN_ADO.PHENO1.glm.linear MAIN_ADO.Ps; mv MAIN_ADOH.PHENO1.glm.linear MAIN_ADOH.Ps; mv MAIN.PHENO1.glm.logistic MAIN.Ps; mv MAIN_AD.PHENO1.glm.logistic MAIN_AD.Ps; mv MAIN_ADO.PHENO1.glm.logistic MAIN_ADO.Ps; mv MAIN_ADOH.PHENO1.glm.logistic MAIN_ADOH.Ps; awk '{print $3,$9,$12}' MAIN.Ps > MAINPs0; awk '{print $3,$9,$12}' MAIN_AD.Ps > MAIN_ADPs0; awk '{print $3,$9,$12}' MAIN_ADO.Ps > MAIN_ADOPs0; awk '{print $3,$9,$12}' MAIN_ADOH.Ps > MAIN_ADOHPs0; sed -i 's/NA/0.99999999999999/' MAINPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOHPs0; sed -i 's/NA/0.99999999999999/' MAINPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOPs0; sed -i 's/NA/0.99999999999999/' MAIN_ADOHPs0; sed '1d' MAINPs0 > MainPs_Add1; sed '1d' MAIN_ADPs0 > MainPs_AddDom1; sed '1d' MAIN_ADOPs0 > MainPs_AddDomOD1; sed '1d' MAIN_ADOHPs0 > MainPs_AddDomODHet1; sed -e 's/ /\t/g' MainPs_Add1 > MainPs_Add; sed -e 's/ /\t/g' MainPs_AddDom1 > MainPs_AddDom; sed -e 's/ /\t/g' MainPs_AddDomOD1 > MainPs_AddDomOD; sed -e 's/ /\t/g' MainPs_AddDomODHet1 > MainPs_AddDomODHet"
            
            ret = subprocess.run(command, capture_output=True, shell=True)
            # before Python 3.7:
            # ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
            print(ret.stdout.decode())
            





print(' "SAFE-h2 has now finalized the combinatory VCF file and corresponding bfiles." ')







print('Cleaning!')




command = "rm MAIN_AD0.bim; rm MAIN_ADO0.bim; rm MAIN_ADOH0.bim; rm head.vcf; rm ADOH_MergePs1; rm ADOH_MergePs2; rm ADOH_MergePs3; rm Heter*; rm MAIN_ADOHPs0; rm MAIN_ADOH.Ps; rm MainPs_AddDomODHet1; rm dominance*; rm overdominance*; rm MainPs_AddDom1; rm MainPs_AddDomOD1; rm MainPs_Add1; rm MAIN.Ps; rm MAIN_AD.Ps; rm MAIN_ADO.Ps; rm MAINPs0; rm MAIN_ADPs0; rm MAIN_ADOPs0; rm *.vcf; rm Additive*; rm Dominance*; rm Overdominance*; rm Add_OD_Dom.txt; rm Prog7_select_act.txt; rm Prog1_select_act.txt; rm Prog2_select_act.txt; rm Prog3_select_act.txt; rm Prog4_select_act.txt; rm Dominance.vcf; rm Overdominance.vcf; rm *.nosex; rm *.tped; rm *.tfam; rm *.map; rm *.BN.kinf; rm *.log; rm Dominance.reml; m Overdominance.reml; rm Additive.reml; rm AD_MergePs1; rm AD_MergePs2; rm AD_MergePs3; rm ADO_MergePs1; rm ADO_MergePs2; rm ADO_MergePs3; rm Additive.vcf; rm AdditiveEx.vcf; rm AdditiveEx1.vcf; rm DominanceEx.vcf; rm DominanceEx1.vcf; rm OverdominanceEx.vcf; rm OverdominanceEx1.vcf; rm Merged_Add-Dom1.vcf; rm Merged_Add-Dom2.vcf; rm Merged_Add-Dom3.vcf; rm Merged_Add-Dom_OD1.vcf; rm Merged_Add-Dom_OD2.vcf; rm Merged_Add-Dom_OD3.vcf; rm headers"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())


command = "mkdir SAFE-h2_output_files; mv Main* SAFE-h2_output_files; mv MAIN_AD* SAFE-h2_output_files; mv MAIN_ADO* SAFE-h2_output_files; mv MAIN_ADOH* SAFE-h2_output_files; mv *Merge* SAFE-h2_output_files"
ret = subprocess.run(command, capture_output=True, shell=True)
# before Python 3.7:
# ret = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
print(ret.stdout.decode())


        
final = pyfiglet.figlet_format('SAFE-h2', width=100)
print(f'[white]{final}[/white]')

print('"has completed the analyses. The SAFE-h2 version is 2026a."')

print('SAFE-h2 Citation: Behrooz Darbani, Mogens Nicolaisen. On the genetic origions of phenotypes in genome-wide association studies: the SAFE-h2 tool for exploring additive and non-additive allelic effects. BMC Bioinformatics 2026. https://doi.org/10.1186/s12859-026-06464-6')

exit()
