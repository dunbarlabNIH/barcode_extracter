import pandas as pd
import numpy as np
import glob
import cPickle
import multiprocessing as mp
import sys
from functools import partial

'''
    Code created by Diego A. Espinoza, Samson J. Koelle, and Rong Lu.
'''


def extractBarcodes(fastq, libID, threshold):
    '''
    this function will create a dictionary
    of the barcodes in a fastq file with
    reads over 100.
    '''
    print "Extracting " + str(fastq)
    file_stream = open(fastq,'r')
    nonthresh_dict = {}
    id_length = len(libID)
    reads = 0
    mapped = 0
    for i, line in enumerate(file_stream):
        if (i % 4 == 1):
            reads += 1
            if line[:id_length] == libID:
                mapped += 1
                if line[:50] in nonthresh_dict:
                    nonthresh_dict[line[:50]] += 1
                else:
                    nonthresh_dict[line[:50]] = 1
    file_stream.close()
    readmeinfo = (fastq, mapped, reads, (100*float(mapped)/reads),  threshold, libID)
    thresh_dict = {}
    for key in nonthresh_dict:
        if nonthresh_dict[key] > threshold:
            thresh_dict[key] = nonthresh_dict[key]
    return (thresh_dict, fastq, readmeinfo)



if __name__ == "__main__":
    '''
    First argument should be either Y or N
    Second argument should be name of file to use if first argument is Y, should be X if first argument is N.
    Third argument should be 6 base-pair library ID for barcode library
    Fourth argument should be length of barcodes in barcode library
    Fifth argument should be sequence read threshold that barcode must pass to be stored for each FASTQ file
    Sixth argument should be desired file name for pickled output of extract_barcodes.py
    Example of usage within bash:
    python extract_barcodes.py N X GTAGCC 35 100 ZH33_extracted.bin
    '''
    
    previous_boolean = sys.argv[1]
    old_file_name = sys.argv[2]
    userlibID = sys.argv[3]
    barcode_length = int(sys.argv[4])
    input_thresh = int(sys.argv[5])
    newfile = sys.argv[6]
    
    
    if(previous_boolean == "Y"):
        old_file = open(old_file_name,'r')
        old_dictlist = cPickle.load(old_file)
        previous_files = [x[1] for x in old_dictlist]
        fileList = []
        for file in glob.glob("*.fastq"):
            if file not in previous_files:
                fileList.append(file)
    elif(previous_boolean == "N"):
        fileList = []
        for file in glob.glob("*.fastq"):
                fileList.append(file)
    else:
        exit("Must choose Y or N.")

    dictpool = mp.Pool(processes=4)
    partial_extract = partial(extractBarcodes, libID = userlibID, threshold = input_thresh)
    new_dictlist = dictpool.map(partial_extract, fileList)
    if(previous_boolean == "Y"):
        new_dictlist = new_dictlist + old_dictlist

    file_dump = open(newfile,'w')
    cPickle.dump(new_dictlist, file_dump)
    file_dump.close()


