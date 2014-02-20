#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 10:13:48 2013

@author: amnon amir

Look at all the OTUs for a given taxa (i.e. Gammaproteobacteria) and select 
candidate contamination sequence OTUs.
output: a list of sorted taxa OTUids according to mean frequency
"""

import sys
import argparse
import numpy as np
from biom.parse import parse_biom_table
from biom.util import biom_open

def get_high_freq_otus(biom_table,taxonomyclass,taxonomyname,level):
    """Return a list of all OTUs from the taxonomy with cumulative freq>level 

    Arguments:
        - biom_table: the biom table to get OTU frequencies from
        - taxonomyclass:    integer describing  the class of the taxonomy to 
                            process (i.e. 0-kingdom,1-phylum etc...)
        - taxonomyname: name of the taxonomy to process 
                        (i.e. 'Gammaproteobacteria)
        - level:    the cumulative frequency above which to return the OTUs
                    (i.e. 0.03 for 3%). If 0, all OTUs in the 
                    taxonomy are sorted and output.
                
    Note:   Output is obtained by filtering only the OTUs which are inside the
            taxonomy requested, then sorting them according to their frequency
            at ascending order, and then calculating the cumulative frequency.
            Once the cumulative freq. reaches 'level', all OTUs beyond this
            are returned (i.e. all the highest freq. OTUs).
            Output format:
            a list of strings with:
            OTUid   \t  MeanFreq    \t  CumulativeFreq
    """

    # normalize the freqs. of each sample
    tablenorm = biom_table.normObservationBySample()

    # keep only otus which are from the desired taxonomy and their mean
    # frequency in all samples
    odat=[]
    class_idx = taxonomyclass
    for values, ids, metadata in tablenorm.iterObservations():
        tname=metadata['taxonomy'][class_idx].lstrip()
        if tname == taxonomyname:
            mv = np.mean(values)
            odat.append((ids,mv))
    
    assert len(odat)>0,'no OTUs matching taxonomy %r at level %d found' %(taxonomyname,class_idx)
    # sort filtered otus according to mean freq. (small to big)
    odat.sort(key=lambda tup: tup[1])

    # calculate the cumulative mean frequency
    csum=[(odat[0][0],odat[0][1],odat[0][1])]
    for cval in odat[1:]:
        csum.append((cval[0],cval[1],csum[-1][2]+cval[1]))
    
    # now get it from big to small
    csum.reverse()
    
    # and return everything above the threshold (for filtering)
    # as a list of tab delimited strings
    result=[]
    for cval in csum:
        if cval[2]>=level:
            result.append(cval[0]+"\t"+str(cval[1])+"\t"+str(cval[2]))
    return(result)


def main(argv):
    parser=argparse.ArgumentParser(description=
        'Select Gammaproteobacteria (or other group) contamination candidates')
    parser.add_argument('-i','--biom',help='biom file of the experiment')
    parser.add_argument('-o','--output',help='output file name')
    parser.add_argument('-c','--classpos',
                        help='class of taxonomy name (0-kingdom,1-phylum etc.',
                        default=2,type=int)
    parser.add_argument('-t','--taxonomy',
                        help='taxonomy name (including c__ or equivalent)',
                        default='c__Gammaproteobacteria')
    parser.add_argument('-l','--level',help=
                        'minimal cumulative level to filter (0 to get all)',
                        default='0.03',type=float)
    
    args=parser.parse_args(argv)

    # load the biom table
    biom_table = parse_biom_table(biom_open(args.biom,'U'))
    # find the high freq. OTUs
    result=get_high_freq_otus(biom_table,args.classpos,args.taxonomy,args.level)
    
    # and write them to the file
    with open(args.output,'w') as snames:
        for cstr in result:
            snames.write(cstr+'\n')

if __name__ == "__main__":
    main(sys.argv[1:])          
