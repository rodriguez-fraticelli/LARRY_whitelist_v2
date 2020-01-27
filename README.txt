This script is used to compare two PCR replicates of the LARRY EGFP v1 barcode library prep.  

The usage is:

python whitelist_script_v2.py File1.fastq File2.fastq Hamming-distance Minimum-reads 

The script takes four arguments (separated by spaces): 
1)	First replicate fastq
2)	Second replicate fastq
3)	Hamming-distance (we use “4” for this size of barcode)
4)	Minimum-reads (we use 5, for 5 million reads per replicate)

The output of the script is a filtered set of fastqs (which grabs only reads with correct barcode structure), two lists of min-reads filtered barcodes, for each first and second replicate files, and one table with the barcodes detected in both replicates and collapsed by hamming distance (whitelist_merged.csv). 

This script is still in development. 
