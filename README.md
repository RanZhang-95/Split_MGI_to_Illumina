# Split_MGI_to_Illumina


This repository contains two Python scripts that allow fastqs read1 and read2 data sequenced by MGI to be splitted based on sample indexes into individual fastq files. These two scripts also convert MGI header into Illumina header so that downstream pipeline Cellranger can be run.

Sample indexes should be provided as a text file. read2 file contains information for sample indexes. Thus split_MGI_read2_to_illumina.py should be run before running split_MGI_read1_to_illumina.py.


split_MGI_read2_to_illumina.py takes a sample indexes text file and read2 fastq.gz file, convert the  MGI header of the fastq file and split records into corresponding sample files. Noted that the sample barcodes are located at the beggining of the reverse complement read2 sequence. The program returns n sample fastq files with changed header and sample indexes attaching to every record ID. It also returns a fastq file containing records that have no sample indexes found. Addiotionally, the program returns a record_barcode.txt to record illumina header. 

split_MGI_read1_to_illumina.py splits read1 fastqs files based on record_barcode.txt.
