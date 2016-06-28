"""
find_diffs.py

Finds sites with an alternate base frequency of at least min_freq.
Frequency is calculated across all samples at a given site.
"""
import sys
from mpileLine import mpileLine
from mpSample import mpSample
import traceback
import subprocess
from pyBinom import pbinom

def min_alt_freq_all_samples(mutant,allSamples,min_freq):
    """

    """
    totalDepth = 0
    otherAltCount = 0
    mutantAltBase = mutant.majorAltBase

    for sample in allSamples:
        if sample != mutant:
            totalDepth += int(sample.depth)
            otherAltCount = otherAltCount + sample.getBaseCount(mutantAltBase[0]) + \
                    sample.getBaseCount(mutantAltBase[1])

    if totalDepth != 0:
        freqAltOther = float(otherAltCount)/float(totalDepth)
    else:
        freqAltOther = 0 

    return freqAltOther > min_freq


def freq_all_samples(mutant_base, all_samples):
    """
    """
    total_depth = 0
    alt_count = 0
    
    for sample in all_samples:
        total_depth += int(sample.depth)
        alt_count += sample.getBaseCount(mutant_base[0]) + \
            sample.getBaseCount(mutant_base[1])

    if total_depth != 0:
        freq = float(alt_count)/float(total_depth)
    else:
        freq = 0

    return freq


if __name__ == "__main__":
    mpile_in,mpile_out = sys.argv[1],sys.argv[2]

    outFile = open(mpile_out,"w")   

#    min_alt = int(sys.argv[3])
    min_freq = 0.3
 
    with open(mpile_in) as mp_file:
        for line in mp_file:
            try:
                mpLine = mpileLine(line)
                if mpLine.chrom not in "pseudo0mitochondrionchloroplast":
                    mutant = mpLine.getMutant()

                    if mutant != None:
                        alt_freq = freq_all_samples(mutant.majorAltBase, mpLine.samples)

                        if alt_freq >= min_freq:
                            outFile.write(line)

                        else:
                            print alt_freq

            except:
                print "ERROR"
                print sys.exc_info()[0] 
                print traceback.format_exc()
                print line
             
    outFile.close()



