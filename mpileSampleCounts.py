import sys
from mpileLine import mpileLine
from mpSample import mpSample
import traceback
import subprocess

new_CC_smpls = ["CC_A","CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O","CC_P"]
new_GP_smpls = ["GP_A","GP_B","GP_C","GP_D","GP_E","GP_F","GP_G","GP_H","GP_I","GP_J","GP_K","GP_L","GP_M","GP_N","GP_O","GP_P"]

old_CC_smpls = ["CC_B","CC_C","CC_D","CC_E","CC_F","CC_H","CC_I","CC_J","CC_K","CC_L","CC_N","CC_O"]
old_GP_smpls = ["GP_B","GP_C","GP_D","GP_E","GP_F","GP_G","GP_H","GP_I","GP_J","GP_K","GP_L","GP_M","GP_N","GP_O"]



idList = None
counts = {}

if __name__ == "__main__":
    mpile_in,out_file,old_new_sample = sys.argv[1],sys.argv[2],sys.argv[3]
    outFile = open(out_file,"w")   
    
    if old_new_sample == "old":
        CC_smpls = new_CC_smpls
        GP_smpls = new_GP_smpls

    if old_new_sample == "new":
        CC_smpls = old_CC_smpls
        GP_smpls = old_GP_smpls
 
    #Assign idList to one of CC or GP smples
    if "CC" in mpile_in:
       idList = CC_smpls 
    elif "GP" in mpile_in:
        idList = GP_smpls
    else:
        print "Error incorrect mpile_in name"
        sys.exit()
    
    ##Loop though mpile file 
    with open(mpile_in) as mp_file:
        for line in mp_file:
            try:
                mpLine = mpileLine(line)
                mutant = mpLine.getMutant()
                mutantID = mpLine.getMutantID(idList) 

                counts[mutantID] = counts.get(mutantID,0) + 1
        
            except:
                print sys.exc_info()[0] 
                print traceback.format_exc()
                print line
    
    #Write out counts dict to outFile          
    for key in counts.keys():
        outFile.write(str(key) + ":" + str(counts[key]) + "\n")
    
         
    outFile.close()



