from mpSample import mpSample
import re

class mpileLine:

    def __init__(self,raw_line):
        self.sline = raw_line.replace(" ","\t").split("\t") 
        self.raw_line = raw_line
        self.chrom = self.sline[0]
        self.pos = self.sline[1]
        self.refBase = self.sline[2]
        self.siteID = self.chrom + ":" + self.pos
        self.mutantSample = None
        self.mutantIDX = None
        self.samples = []
        self.loadSamples()
        self.loadMutant()
        return
    
    def getMutantID(self, IDlist):
        """
        Returns the sample name of the mutant given a list of 
            sample names corresponding to the order ran through mpileup.
        """
        if self.mutantIDX == None:
            return "None"

        return IDlist[self.mutantIDX]


    def loadSamples(self):
        i = 3 
        while i < len(self.sline):
            sample = mpSample(self.sline[i:i+4])
            self.samples.append(sample)
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
#        max_reads = 0
#        for sample in self.samples:
#
#            if sample.majorAltBaseCount > max_reads:
#                print str(sample.majorAltBaseCount) + " > " + str(max_reads)
#                self.mutantSample = sample
#                max_reads = sample.majorAltBaseCount
#            else:
#                print str(sample.majorAltBaseCount) + " < " + str(max_reads)

        max_reads = 0
        for i in range(0,len(self.samples)):
            if self.samples[i].majorAltBaseCount > max_reads:
                max_reads = self.samples[i].majorAltBaseCount
                self.mutantSample = self.samples[i]
                self.mutantIDX = i
        
        #if self.mutantSample == None:
        #    print "No alternate reads found in any sample"        

        return


    def repr(self):
        """
        Returns a string representation of the mpileLine, including the samples
        """
        repr_str = ""
        
        repr_str = repr_str + self.raw_line + "\n"
        for s in self.samples:
            repr_str = repr_str + s.repr()

        return repr_str 



