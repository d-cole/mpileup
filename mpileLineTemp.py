from mpSample import mpSample
import re

class mpileLineTemp:

    def __init__(self,raw_line):
        self.sline = raw_line.replace(" ","\t").split("\t") 
        self.raw_line = raw_line
        self.chrom = self.sline[0]
        self.pos = self.sline[1]
        self.refBase = self.sline[2]
        self.siteID = self.chrom + ":" + self.pos
        self.mutantSample = None
        self.samples = []
        self.loadSamples()
        self.loadMutant()
        return

    def loadSamples(self):
        i = 3 
        while i < len(self.sline):
            sample = mpSample(self.sline[i:i+4])
            self.samples.append(sample)
#            print sample.__repr__()
            i += 4
        return 

    def getSample(self,pos):# -> mpileSample
        """
        Returns the sample at the given pos
        Used to coordinate samples between vcf and pileup files
        Assumes samples arranged in the same order
        """
        return self.samples[pos]

    def getMutant(self):
        """
        Returns the mpSample object for the 'mutant'
        """
        return self.mutantSample

    def loadMutant(self):
        """
        Assigns self.mutantSample to the sample with the most alt reads
        """
        max_reads = 0
        for sample in self.samples:

            if sample.majorAltBaseCount > max_reads:
#                print str(sample.majorAltBaseCount) + " > " + str(max_reads)
                self.mutantSample = sample
                max_reads = sample.majorAltBaseCount
#            else:
#                print str(sample.majorAltBaseCount) + " < " + str(max_reads)

    def repr(self):
        """
        Returns a string representation of the mpileLine, including the samples
        """
        repr_str = ""
        
        repr_str = repr_str + self.raw_line + "\n"
        for s in self.samples:
            repr_str = repr_str + s.repr()



        return repr_str 



