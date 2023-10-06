Pipeline outline: 
1. Build your environments. 
2. Get your barcode file. 
3. Run the analysis. 
3.5. QC mapping percentage with mapping stats output from samtools. 
3.75. QC demix coverage (make sure it is above 80%)
4. Frankenstein your dashboard files using the demix file and asset tiger in excel, save as tab deliminated txt file. 
5. Remove windows line endings. 
6. Upload data onto data depot. 



Other notes:
To run the covid_analysis.sh pipeline, need to build three environments: 
1) the freyja environment, 
2) the fastqc environment for qc purposes. Fastqc doesnt play well with freyja, so they need to be separate and, 
3) the pre-processing environment used for trimming. 

"Freyja update" is not working appropriately right now. 
In order to use the most up to date freyja barcode, go to the freyja data repository and download the .csv barcode file of you choice under history_barcodes. 
Change the variable in the script to reflect the location of this barcode file. 

*** The script runs each barcode individually. You will need to submit a script for each barcode ***

Update the variables at the beginning of the script, and the paths to the conda environments (at lines 35, 45, 51, 57, 78).

Output files will be written to a folder you specify in the variable $WORKINGDIR

Three .yaml files for building conda environments:


For freyja:
name: freyja_env_scratch
channels:
  - defaults
  - bioconda
  - conda-forge
dependencies:
  - python=3.7
  - freyja
  - ivar
  - samtools
  - UShER
  - cvxpy
  - numpy
  - pandas

For fastqc: 
name: fastqc_env_scratch
channels:
  - bioconda
dependencies:
  - fastqc
  
For pre-processing, etc: 
name: pre_process_env
channels:
  - defaults
  - bioconda
dependencies: 
  - python=3.7
  - pycoqc
  - nanoplot
  - blast
  - minimap2
  - trimmomatic
  - samtools
  


