"""
find_excl_other_genotype.py

Given a .mpileup file of sites, find these sites in other_mpileup 
    and write them to a file.
"""
import sys
from mpileLine import mpileLine
from mpSample import mpSample

if __name__ == "__main__":

    focal_excl_sites = sys.argv[1]
    other_mpileup = sys.argv[2]

    out_file = open(sys.argv[3], "w")

    #Load in excl sites from focal genotype
    # Presence of a site indicated by 1 at sites[chrom][pos]
    # sites format:
    #   {CHROM:{POS:INT, POS:INT, ...}, CHROM:{...}, ...}
    sites = {}

    with open(focal_excl_sites) as f:
        for line in f:
            sline = line.split()
            sites[sline[0]] = sites.get(sline[0],{})
            sites[sline[0]][sline[1]] = 1
 
    
    #Slow part
    #Check every line in other mpileup file looking for this site
    with open(other_mpileup) as fmp:
        for line in fmp:
            sline = line.split()

            #Check if chrom matches
            if sites.get(sline[0], None) != None:

                #Check if pos present in sites
                if sites[sline[0]].get(sline[1], 0) == 1:
                    out_file.write(line)


    out_file.close()




