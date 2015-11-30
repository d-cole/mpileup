import sys
from mpileLineTemp import mpileLineTemp
from mpSample import mpSample
import traceback
import subprocess

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


    binomResult = "TRUE" in str((subprocess.Popen("~/spirodela/r/binom.r " + str(altReads)\
             +" "+str(refReads),shell=True,stdout=subprocess.PIPE)).communicate()[0])

    return binomResult


if __name__ == "__main__":
    mpile_in,mpile_out = sys.argv[1],sys.argv[2]
    outFile = open(mpile_out,"w")   
 
    with open(mpile_in) as mp_file:
        for line in mp_file:
            try:
                mpLine = mpileLineTemp(line)
                mutant = mpLine.getMutant()
                
                if mpLine.chrom not in "pseudo0mitochondrionchloroplast":
    #                print mutant.__repr__()
                        if mutant.altReadCount(30) >= 3:
                            if mutant.altReadBothStrands(30):
                            #if mutant.altReadBothStrands():
                                if maxFreqAltOtherSamples(mutant,mpLine.samples,0.02):
                                    #print "pass max freq other samples"
                                    if passBinom(mutant):
                                        #print "passed binom test
                                        outFile.write(line)
        
            except:
                print sys.exc_info()[0] 
                print traceback.format_exc()
                print line
             
            
    outFile.close()



