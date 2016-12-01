#barcode_extracter

This is python code used to extract cellular barcoding data from FASTQ files.
Final output is a table with cellular barcode sequences as row names, and FASTQ files as columns and counts corresponding to each barcode present within that FASTQ file.

Code by Diego A. Espinoza, Samson J. Koelle, and Rong Lu.
Created at Laboratory of Cynthia E. Dunbar at National Heart, Lung, and Blood Institute at National Institutes of Health, Bethesda, MD.

Please use the following steps to extract cellular barcode data from FASTQ files:

1. Place all FASTQ files with data to extract in a single directory.
2. Place "extract_barcodes.py" in the same directory.
3. Run "extract_barcodes.py" within that directory as shown:
    a. NOTE: If a *.bin file exists with extracted FASTQ file data, you can add new data to this existing extraction. To do this, the existing *.bin file should be placed in same directory as FASTQ files and code, and BOOLEAN and OLDFILE options should be used correctly as specified below. Adding data to existing extraction in this way excludes re-extraction of any FASTQ files in existing *.bin extraction that may be in the directory.
```
python extract_barcodes.py BOOLEAN OLDFILE LIBID BARCODELENGTH THRESH OUTFILENAME
BOOLEAN             = Y if existing *.bin file with extracted data exists and need to add new data to old *.bin file
                      N if no *.bin file with extracted data exists
OLDFILE             = if BOOLEAN is Y, must be name of *.bin file in directory with extracted data
                      if BOOLEAN is N, use X as argument
LIBID               = 6 base pair DNA sequence of correct library ID
BARCODELENGTH       = length of barcodes in library
THRESH              = sequence read threshold that barcode must pass to be stored for each FASTQ file
OUTFILENAME         = desired file name for pickled output of extract_barcodes.py

Example usage 1:
python extract_barcodes.py N X GTAGCC 35 100 ZH33_extracted.bin

Example usage 2:
python extract_barcodes.py Y ZH33_extracted_20160730.bin GTAGCC 35 100 ZH33_extracted_20161130.bin
```

4. Once extract_barcodes.py is finished, run "combine_barcodes.py" in same directory where pickled output of "extracted_barcodes.py" is located.
    a. NOTE: In order to run efficiently on cluster systems, exclusive node allocation is necessary to preform multi-core processing.

```
python combine_barcodes.py EXTRACTEDFILENAME FINALFILENAME LIBID BARCODELENGTH MAXMISMATCHES CORES
EXTRACTEDFILENAME   = file name of extract_barcodes.py output
FINALFILENAME       = desired file name for text output of combine_barcodes.py
LIBID               = 6 base pair DNA sequence of correct library ID
BARCODELENGTH       = length of barcodes in library
MAXMISMATCHES       = maximum number of indels/mismatches to use for combining barcodes
CORES               = number of cores to use for processing

Example usage:
python combine_barcodes.py ZH33_extracted.bin ZH33_combined.txt GTAGCC 35 2 30
```


5. Upon completion, output of "combine_barcodes.py" will be a table that can be used for analysis.