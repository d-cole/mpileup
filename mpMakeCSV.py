import sys
from mpileLine import mpileLine
from mpSample import mpSample
import traceback
import subprocess


if __name__ == "__main__":
    mpile_in,csv_out = sys.argv[1],sys.argv[2]
    outFile = open(csv_out,"w")   
    
    columns = "pos,refCount,altCount"
    outFile.write(columns + "\n")
 
    with open(mpile_in) as mp_file:
        for line in mp_file:
            try:
                mpLine = mpileLine(line)
                mutant = mpLine.getMutant()
                outString = mpLine.chrom + ":" + mpLine.pos + ","
                outString = outString + str(mutant.refBaseCount) +"," + str(mutant.majorAltBaseCount)
                outFile.write(outString + "\n")
        
            except:
                print sys.exc_info()[0] 
                print traceback.format_exc()
                print line
             
            
    outFile.close()



