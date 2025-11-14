#################################
SEQSARS
#################################

The subtf_seqsars.sh script is the most common wrapper script for the SEQSARS workflow, as 
it is most often run on Thorny Flat. The workflow accepts sequence read data as input 
(fastq.gz format) and generates COVID variant proportions as output. This output can be piped 
directly into the SEQR workflow to generate files for the state dashboard and NCBI repositories.



*** SETUP ***

* First, make sure that you have configured all of the required conda environments! Config 
files (yml format) are provided in the resources subdir to assist with this process.

* Log into the Thorny Flat HPC, and submit a job with the subtf_seqsars.sh script to execute. 



*** PARAMETERIZATION ***

The one *required* argument to subtf_seqsars.sh is:

(-i) Path to the run directory. This is the parent run folder that contains the 'fastq_pass' 
folder (among others).


There are two additional parameters that we *strongly recommended* you also provide:

(-r) A run ID. This is usally the very top-level folder name (eg, SARS_20251105). By default, 
I will create a runid for you that includes the date of the script execution. This is less 
than good and is discouraged.

(-o) Path to a valid output directory. This folder will contain all of the output files, 
including copies of the read files, so it will require some amount of space. By default, I 
will create folder called 'output/[RUNID]' in the run directory. This is *strongly discouraged*. 


There are two more parameters that you can customize, but rarely need to do so:

(-g) Path to a reference genome file. The default file is located in the project resources 
subfolder and is named 'NC_045512_Hu-1.fasta'.

(-b) Path to an usher barcodes file. The default file is located in the project resources 
subfolder and is named 'usher_barcodes.csv'. This file contains the barcodes (SNP combinations) 
used by freyja to assign variant ids. If it is missing, you can try to download it from here:

	https://github.com/andersen-lab/Freyja/blob/main/freyja/data/usher_barcodes.csv



*** EXECUTION ***

The subtf_seqsars.sh script initializes the conda env variables and calls the sh/seqsars.sh 
script. If you aren't on Thorny Flat, you can simply init the conda env variables and run 
sh/seqsars.sh directly - it accepts all the same parameters as the subtf_seqsars.sh script.

The sh/seqsars.sh script is the primary controller for the SEQSARS workflow. It creates the 
output directory structure and calls out to the required packages (minimap2, freyja, etc.) 
to do the heavy lifting of the workflow.

At the end of the SEQSARS workflow, the output directory will contain the following subdirs:
	0 RUN_FILES
	1 BASECALLING
	2 FQ_CONCAT
	3 QC1_RESULTS
	4 FQ_CHECKED
	5 QC2_RESULTS
	6 MAPPED
	7 DEMIX
	8 DASHBOARD
	9 SRA

