import sys
from mpileLineTemp import mpileLineTemp
from mpSample import mpSample
#
#TE = {Pseudo0:[(x,y),(z,l)...],Pseudo1:[(m,x)...]...}
TE_ranges = {}
CHROM,START,STOP = 0,1,2

def loadTEranges(TE_file_loc):
    """
    Load TE ranges from a text file of the format CHROM START STOP
    """
    with open(TE_file_loc) as TE_file:
        for line in TE_file:
            line_col = str.split(line)
            #Identify chromosome
            TE_ranges.setdefault(line_col[CHROM],[]).append((line_col[START],line_col[STOP]))
    TE_file.close()
    return

def validRange(chrom,pos):
    """
    Determine if the given site is withen any TE range
    """
# any(lower <= postcode <= upper for (lower, upper) in [(1000, 2249), (2555, 2574), ...])
    if any(float(low) <= float(pos) <= float(high) for (low,high) in TE_ranges[chrom]):
        return False
    return True

if __name__ == "__main__":
    TE_file_loc,mp_file_loc = sys.argv[1],sys.argv[2]
    loadTEranges(TE_file_loc)

    trimmed_mp = open(mp_file_loc[:mp_file_loc.find(".")] + "_TEremoved.mpileup","w")
    removed_sites = open(mp_file_loc[:mp_file_loc.find(".")] + "_TE_sites.mpileup","w")

    with open(mp_file_loc) as mp_file:
        for line in mp_file:
            
            mp_line = mpileLineTemp(line)
            if validRange(mp_line.chrom,mp_line.pos):
                trimmed_mp.write(line)
            else:
                trimmed_mp.write(line)


    removed_sites.close()
    mp_file.close()
    trimmed_mp.close()
