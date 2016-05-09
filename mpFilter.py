import sys
from mpileLine import mpileLine
from mpSample import mpSample
import traceback
import subprocess
from pyBinom import pbinom

def maxFreqAltOtherSamples(mutant,allSamples,freq):
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

#    print "Depth: " + str(totalDepth)   
#    print "AltCount: " + str(otherAltCount)
#    print "freq: " + str(freqAltOther)
#
    return freqAltOther < freq


def passBinom(mutant):
    """
    
    """
    altReads = mutant.majorAltBaseCount
    refReads = mutant.refBaseCount
#    print "alt: " + str(altReads)
#    print "ref: " + str(refReads)

    return pbinom(altReads,refReads,0.5) > 0.02
#    binomResult = "TRUE" in str((subprocess.Popen("~/spirodela/r/binom.r " + str(altReads)\
#             +" "+str(refReads),shell=True,stdout=subprocess.PIPE)).communicate()[0])
#
#    return binomResult
#

if __name__ == "__main__":
    mpile_in,mpile_out = sys.argv[1],sys.argv[2]
#    no_mutant = sys.argv[3]

    outFile = open(mpile_out,"w")   
    #no_mut = open(no_mutant,"w")    

    min_alt = int(sys.argv[3])
 
    with open(mpile_in) as mp_file:
        for line in mp_file:
            try:
                mpLine = mpileLine(line)
                mutant = mpLine.getMutant()
                
                if mpLine.chrom not in "pseudo0mitochondrionchloroplast":

                        #No mutant at this site
                       # if mutant == None:
                       #     no_mut.write(line)
                       #     continue

                        if mutant.altReadCount(30) >= min_alt:
                            if mutant.altReadBothStrands(30):

                                if maxFreqAltOtherSamples(mutant,mpLine.samples,0.02):

                                    if passBinom(mutant):

                                        outFile.write(line)
        
            except:
                print "ERROR"
                print sys.exc_info()[0] 
                print traceback.format_exc()
                print line
             
    #no_mut.close() 
    outFile.close()



