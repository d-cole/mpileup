"""
mpileSample.py 
Class for representing samples from the .pileup format
"""

base_options = ["Aa","Tt","Gg","Cc"]


class mpileSample:
    """
    ATCGN - alt base on forward strand
    atcgn - alt base on reverse strand
    .   - ref on forawrd strand
    ,   - ref on reverse strand
    \+[0-9]+[ACGTNacgtn]+   - insertion of one or more bases
    -[0-9]+[ACGTNacgtn]+    - deletion of one or more bases
    ^   - start of read segment ACII char after '^' - 33 == mapping quality
    $   - end of a read segment
    *   - placeholder for del base in multiple bp deletion
    <   - reference skip
    >   - reference skip
    """

    def __init__(self,sampleSegment):
        self.readCount = sampleSegment[0]
        self.baseString = sampleSegment[1]
        self.baseQuality = sampleSegment[2]
        self.mappingQual = sampleSegment[3]#Not sure about this
        self.baseString_masked = ""
        
        self.baseCounts = {"A":[],"a":[],"T":[],"t":[],"C":[],"c":[],"G":[],"g":[],"N":[],"n":[],".":[],",":[]}
        self.majorAltCount = 0 
        self.majorAltBase = "NA"
        self.refCount = 0

        self.maskIndels()
        self.parseBaseString()
        self.getMajorAlt()
    #======================= NOT TESTED YET ================ 
    def altReadExist(self,min_alt,min_qual = 0):
        """
        Returns where there are >= min_alt alternate reads
        ADD MIN REQ ON BASE QUALITIES
        """
        
        if self.majorAltBase in base_options:
            good_qual_count = 0
            #baseCounts stores the indexes of the base occurances
            #the majorAltBase string "aA","tT"... is the index to the counts 
            #Loop through all the positions where there is an alt base in the corrresponding base string
            for base_index in (self.baseCounts[self.majorAltBase[0]] + self.baseCounts[self.majorAltBase[1]]):
                 
                #print self.baseString
                #print self.baseQuality
                #print "base_idx: " + str(base_index)

                if (ord(self.baseQuality[base_index]) - 33) >= min_qual:
                    good_qual_count += 1

            return good_qual_count >= min_alt 
        return False
    
    def altReadBothStrands(self,min_qual = 0):
        """
        Returns true if the alt base is present in both directions
        """
        if self.majorAltBase in base_options:
            #get lists of forward and reverse base positions
            forward_bases = self.baseCounts[self.majorAltBase[0]]
            reverse_bases = self.baseCounts[self.majorAltBase[1]]
        
            #Determine if one read exists on each direction
            if len(forward_bases) > 0 and len(reverse_bases) > 0:
                #Determine if there is atleast one high quality read in each direction
                if any((ord(self.baseQuality[base_index]) - 33) >= min_qual for base_index in forward_bases):
                    if any((ord(self.baseQuality[base_index]) - 33) >= min_qual for base_index in reverse_bases):
                        return True

        elif self.majorAltBase != "NA":
            print "mystery base"
            print self.majorAltBase

        return False

    def repr(self):
        """
        Returns the string representation of the mpileSample object
        """
        repr_str = ""
        repr_str = repr_str + "readCount: " + str(self.readCount) + "\n" + \
            "baseString: " + str(self.baseString) + "\n" + \
            "baseQuality: " + str(self.baseQuality) + "\n" + \
            "mappingQual: " + str(self.mappingQual) + "\n" + \
            "MajorAltBase: "  + str(self.majorAltBase) + "\n" + \
            "MajorAltCount: " + str(self.majorAltCount) + "\n" + \
            "refCount: " + str(self.refCount) + "\n\n"

        return repr_str         
    
    def parseBaseString(self):
        """
        Parses the reads from the base string. Stores the information in self.baseCounts
        """
        for i in range(len(self.baseString_masked)):
            char = self.baseString_masked[i]
            if char != "i":
                self.baseCounts[char] = self.baseCounts.get(char,[])
                self.baseCounts[char].append(i)

    def getMajorAlt(self):
        """
        CHANGE METHOD NAME

        Determines the highest frequency ALT allele.
        Sets self.majorAltCount to the number of the major ALT reads.
        Sets self.majorAltBase to the base with the highest read count.
        """
        refCount = len(self.baseCounts["."]) + len(self.baseCounts[","])
        self.refCount = refCount
        aCount = len(self.baseCounts["A"]) + len(self.baseCounts["a"])
        tCount = len(self.baseCounts["T"]) + len(self.baseCounts["t"])
        gCount = len(self.baseCounts["G"]) + len(self.baseCounts["g"])
        cCount = len(self.baseCounts["C"])  + len(self.baseCounts["c"])

        baseCounts = [aCount,tCount,gCount,cCount]
        
        #find max occurance of an alt allele
        countMax = max(baseCounts)
        self.majorAltCount = countMax

        #If max is 0 there are no alternate alleles present
        if countMax != 0:
            for i in range(len(baseCounts)):
                if baseCounts[i] == countMax:
                    self.majorAltBase = base_options[i]
                    break

    def __getIndelLength(self,baseString,indelIndicatorPos):
        """ 
        Given the index to a "+/-" return the length of the indel 
        """
        print "iipos: " + str(indelIndicatorPos)
        print baseString[indelIndicatorPos + 1:] 
        lenString = ""
        for char in baseString[indelIndicatorPos + 1:]:
            if str.isdigit(char):
                lenString = lenString + char
            else:
                break
        print "lenString: " + lenString
        return int(lenString)
        

    def maskIndels(self):
        """
        Creates self.baseString_remIndels from self.baseString.
        mask indel reads with "i"  for easier parsing of reads at this site. 
        Maintains positions of bases for matching with base quality string.
        
        Temporary?
        REMOVE ^x from base segment 
        REMOVE $ from base segment

        """
        print self.majorAltBase
        print self.baseString
        print self.baseString_masked
        print "base count: " + str(self.baseCounts)
        print "baseQuality: " + str(self.baseQuality)        
        print "baseQual len: %s, maskedString len: %s" %(str(len(self.baseQuality)),str(len(self.baseString_masked)))
        

        #Remove insertions & deletions for easier parsing of base string
        i = 0
        while i < len(self.baseString):
            print self.baseString[i]
            if self.baseString[i] == "^":
                #self.baseString_masked = self.baseString_masked + "ii"
                #move over '^x' characters
                i += 2
            else:

                if self.baseString[i] == "+" or self.baseString[i] == "-":
                    #Marks start of indel
                    #indelLength = int(self.baseString[i + 1])
                    indelLength = self.__getIndelLength(self.baseString,i)
                    print "indelLength" + str(indelLength) 
                    #pass over indel 
                    i += indelLength + 2 #add 2 to pass over +/- and onto the read after the indel
#                    self.baseString_masked = self.baseString_masked + "i" * (indelLength + 2)
                else:
                    print "base here: " + self.baseString[i]
                    print "masked: " + self.baseString_masked
                    #Not an indel read
                    if self.baseString[i] != "$":
                        self.baseString_masked = self.baseString_masked + self.baseString[i]
                    i+=1
        return

        



