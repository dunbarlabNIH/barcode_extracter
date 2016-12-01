import pandas as pd
import numpy as np
import glob
import sys
import copy
import cPickle
import multiprocessing as mp
from functools import partial


'''
    Code created by Diego A. Espinoza, Samson J. Koelle, and Rong Lu.
'''


def zR3(seq1, seq2, seq1ind, seq2ind, n_err, totalmatched, barcodelength, maxmismatch):
    '''
    this function takes two sequences (seq1 and seq2) and the start indices (seq1ind and seq2ind)
    and current number of errors (n_err) along with the number of correct bases so far (totalmatched)
    and length of barcode (barcodelength) and tries to match the next index until totalmatched
    matches are found. Uses recursion to return True if two barcodes match with <= maxmismatch
    mismatch/indels, False otherwise
    '''
    if n_err > maxmismatch:
            exit("Critical Error: zR3 saw more than maxmismatches as OK")
    if int(totalmatched) == int(barcodelength):
        return True
    if seq1[seq1ind] == seq2[seq2ind]:
        return zR3(seq1,seq2,seq1ind+1, seq2ind+1,n_err,totalmatched+1, barcodelength, maxmismatch)
    else:
        n_err += 1
        if n_err > maxmismatch:
            return False
        return zR3(seq1, seq2, seq1ind+1, seq2ind, n_err, totalmatched, barcodelength, maxmismatch) or \
               zR3(seq1, seq2, seq1ind, seq2ind+1, n_err, totalmatched, barcodelength, maxmismatch) or \
               zR3(seq1, seq2, seq1ind+1, seq2ind+1, n_err,totalmatched+1, barcodelength, maxmismatch)




def all_disjoint(sets):
    '''
    this function checks to see if the given list of sets are pairwise disjoint
    '''
    all = set()
    for s in sets:
        for x in s:
            if x in all:
                return False
            all.add(x)
    return True





def combine_barcodes(initial_tuple, libID, barcode_length, max_mismatches):
    '''
    takes in a dictionary of barcodes (keys) and their counts (values),
    and uses zr3 to combine the barcodes in the dictionary and returns a combined
    dictionary of barcodes
    '''
    print "Combining barcodes in " + str(initial_tuple[1])
    initial_dict = initial_tuple[0]
    all_barcodes = initial_dict.keys()
    original_barcode_number = len(all_barcodes)
    listofseqs_as_sets = [{i} for i in all_barcodes]
    for i in range(0, len(all_barcodes)):
        for j in range(i + 1, len(all_barcodes)):
            if zR3(all_barcodes[i], all_barcodes[j], len(libID), len(libID), 0, 0, barcode_length, max_mismatches):
                listofseqs_as_sets[i].update([all_barcodes[j]])
    disjoint_set_barcodes = make_sets_disjoint(listofseqs_as_sets, all_barcodes)
    combined_dict = {}
    for setofbar in disjoint_set_barcodes:
        consensus_sequence = max(setofbar, key=lambda i: initial_dict[i])
        subset_dict = dict([(i, initial_dict[i]) for i in setofbar])
        consensus_sum = sum(subset_dict.values())
        combined_dict[consensus_sequence] = consensus_sum
    print "Done combining barcodes in " + str(initial_tuple[1])
    new_barcode_number = len(combined_dict)
    combine_barcodes_info = (initial_tuple[1],) + initial_tuple[2] + (original_barcode_number, new_barcode_number)
    return (combined_dict, initial_tuple[1], combine_barcodes_info)


def make_sets_disjoint(listofsets, listofbarcodes):
    '''
    takes in a list of sets and the list of barcodes and combines
    the list of sets into a set of disjoint sets... derived from
    the theory of connected components in graphs
    '''
    listy = copy.deepcopy(listofsets)
    for barcode in listofbarcodes:
        disjointlistofsets = []
        newset = set()
        for setty in listy:
            if barcode in setty:
                newset = newset.union(setty)
            else:
                disjointlistofsets.append(setty)
        disjointlistofsets.append(newset)
        listy = disjointlistofsets
    assert all_disjoint(listy)
    return listy

def readme_creator(run_info, new_file, libID, barcode_length, max_mismatches):
    '''
    Outputs a README file denoting valuable statistics from the run.
    '''
    readme_file_name = new_file.split(".")[0]+"_README.txt"
    readme_stream = open(readme_file_name, 'w')
    readme_stream.write("Library ID used: " + str(libID) + '\n')
    readme_stream.write("Barcode length used: " + str(barcode_length) + '\n')
    readme_stream.write("Max mismatches used: " + str(max_mismatches) + '\n')
    readme_stream.write("FILENAME" + "\t" + "MAPPED" + "\t" +
                        "READS" + "\t" + "MAPPED %" + "\t" + "Pre-Combined Barcode Number" +
                        "\t" + "Combined Barcode Number" + '\n')
    for my_tuple in run_info:
        readme_stream.write('\t'.join(str(x) for x in my_tuple) + '\n')
    readme_stream.close()



if __name__ == "__main__":
    '''
    First argument should be name of pickled file obtained from extract_barcodes.py
    Second argument should be desired name of output .txt file
    Third argument should be 6 base-pair library ID for barcode library
    Fourth argument should be length of barcodes in barcode library
    Fifth argument should be max number of mismatches to use when combining barcodes
    Sixth argument should be number of cores used for code
    Example of usage within bash:
    python combine_barcodes.py ZH33_extracted.bin ZH33_outfile.txt GTAGCC 35 2 30
    '''
    extracted_file = sys.argv[1]
    user_new_file = sys.argv[2]
    user_libID = sys.argv[3]
    user_barcode_length = int(sys.argv[4])
    user_max_mismatches = int(sys.argv[5])
    user_cores = int(sys.argv[6])
    
    temp = open(extracted_file,'r')
    uncombined_dictlist = cPickle.load(temp)
    temp.close()
    
    second_dictpool = mp.Pool(processes=user_cores)
    partial_combine = partial(combine_barcodes, libID = user_libID, barcode_length = user_barcode_length, max_mismatches = user_max_mismatches)
    combined_dictlist = second_dictpool.map(partial_combine, uncombined_dictlist)
    print "Combining done. Creating new outfile..."
    combined_barcode_dict = [x[0] for x in combined_dictlist]
    combined_barcode_names = [x[1] for x in combined_dictlist]
    combined_barcode_info = [x[2] for x in combined_dictlist]
    new_data_frame = pd.DataFrame.from_records(combined_barcode_dict, index = combined_barcode_names).transpose()
    new_data_frame[np.isnan(new_data_frame)] = 0
    new_data_frame.to_csv(user_new_file, sep = '\t')
    print "Outfile created. Writing info to README..."
    readme_creator(combined_barcode_info, user_new_file, user_libID, user_barcode_length, user_max_mismatches)
    print "README created. Controller done."




