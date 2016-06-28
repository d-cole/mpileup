"""
find_good_diffs.py

Given two .mpileup files (one from each genotype), finds sites where 
    one genotype has an alt base freq > 0.9, the other an alt base freq of 0.

! Would be extremely slow if ran on unfiltered .mpileup
    Further instructions in README
"""
import sys
from mpileLine import mpileLine
from mpSample import mpSample
from find_diffs import freq_all_samples

class site_info:
    """
    Maintains necessary info about a given site.
    A trimmed down mpileLine object
    """
    
    def __init__(self, mpileLine):
        self.chrom = mpileLine.chrom
        self.pos = mpileLine.pos
        self.ref = mpileLine.refBase
        self.avg_depth = self._get_avg_depth(mpileLine.samples)
        self.mutant = mpileLine.getMutant()

        if self.mutant != None:
            self.freq = self._get_freq(mpileLine.samples)
        else:
            self.freq = 0


    def _get_avg_depth(self, samples):
        """
        """
        total = 0.0
        num_samples = 0.0
        for sample in samples:
            total += float(sample.depth)
            num_samples += 1.0

        return total/num_samples


    def _get_freq(self, all_samples):
        """
        Returns frequency of alt base at this site
        """
        total_depth = 0
        alt_count = 0
        
        for sample in all_samples:
            total_depth += int(sample.depth)
            alt_count += sample.getBaseCount(self.mutant.majorAltBase[0]) + \
                sample.getBaseCount(self.mutant.majorAltBase[1])

        if total_depth != 0:
            freq = float(alt_count)/float(total_depth)

        else:
            freq = 0
   
        return freq

    def get_str(self):
        """
        """ 
        site_str = self.chrom + "\t" + self.pos + \
            "\t" + self.ref + "\t" + self.mutant.majorAltBase[0] +\
                "\t" + str(self.freq) +  "\n"

        return site_str 


def load_file_to_dict(file_loc, sites):
    """
    Dictionary format:
    sites: {CHROM: {pos:(CC_info, GP_info), pos: (...)}, CHROM:{...}}
    """
    with open(file_loc) as f:
           for line in f:
               mpile_site = mpileLine(line)
               sites[mpile_site.chrom] = sites.get(mpile_site.chrom, {})
               sites[mpile_site.chrom][mpile_site.pos] = \
                   sites[mpile_site.chrom].get(mpile_site.pos, [])
    
               sites[mpile_site.chrom][mpile_site.pos].append(site_info(mpile_site))
    return 

#! Functionality moved to get_str() in class site_info
#def get_site_str(alt_info):
#    """
#    Returns string representation of a site_info
#    """
#    site_str = ""
#    site_str = site_str + alt_info.chrom + "\t" + alt_info.pos + \
#        "\t" + alt_info.ref + "\t" + alt_info.mutant.majorAltBase[0] +\
#             "\t" + str(alt_info.freq) +  "\n"
#
#    return site_str 

def get_good_sites(sites):
    """
    Given dictionary of sites returns string representation of sites where one genotype has an alt freq > 0.9, 
        the other an alt freq of 0. 

    Requires average depth of both genotypes to be >= 10
    """
    good_sites = []

    for chrom in sites.keys():
        for pos in sites[chrom].keys():
            if len(sites[chrom][pos]) == 2:
                alt_info = sites[chrom][pos][0]
                other_info = sites[chrom][pos][1]

                if alt_info.avg_depth >= 10 and other_info.avg_depth >= 10:
                    if other_info.mutant == None and alt_info.freq > 0.9:
                        good_sites.append(alt_info.get_str())

                
    return good_sites                    
        
    
if __name__ == "__main__":
    alt_excl_loc, other_loc, out_loc = sys.argv[1], sys.argv[2], sys.argv[3]
  
    #sites: {CHROM: {pos:(CC_info, GP_info), pos: (...)}, CHROM:{...}}
    sites = {}
   
    #Load all sites form other genotype file at positions shared with excl_CC 
    #ORDER OF LOADING MATTERS
    # Puts alt_excl genotype sites 
    load_file_to_dict(alt_excl_loc, sites)
    load_file_to_dict(other_loc, sites)

    out_file = open(out_loc, "w")

    good_sites = get_good_sites(sites) 
    for site in good_sites:
        out_file.write(site)

    out_file.close()
        





         
     
