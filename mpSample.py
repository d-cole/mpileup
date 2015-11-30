"""
mpSample.py
Class for representing sample information from the .pileup format
"""
import operator
#base_options stores valid types of bases
base_options = ["Aa","Tt","Gg","Cc"]


class mpSample:
    """
    """

    def __init__(self,sampleSegment):
        """
        sampleSegment is a list of strings containing information pertaining to one sample
        """
        self.bases_masked = ""
        self.majorAltBase = "NA"
        self.majorAltBaseCount = 0
        self.refBaseCount = 0
        self.basePositions = {"A":[],"a":[],"T":[],"t":[],"C":[],"c":[],"G":[],"g":[],"N":[],"n":[],".":[],",":[]}

        self.depth = sampleSegment[0]
        self.bases = sampleSegment[1] 
        self.baseQualities = sampleSegment[2]
        self.mappingQual = sampleSegment[3]

        self.__maskBaseString()
        self.__getBasePositions()
        self.__getMajorAltBase()
    
    def __getIndelLength(self,baseString,indelIndicatorPos):
        """
        Returns the length of the indel
        """
        lenString = ""
        for char in baseString[indelIndicatorPos + 1:]:
            if str.isdigit(char):
                lenString = lenString + char
            else:
                break

        return int(lenString)

    def __maskBaseString(self):
        """
        Creates self.bases_masked from self.bases
        Removes Insertions,deletions,read quality scores.
        len(self.bases_masked) must equal len(self.baseQualities)
        """
        base_pos = 0
        validBases = "aAtTcCgG.,nN*" #use this????
        
        while base_pos < len(self.bases):
            #Order of checks??????
            if self.bases[base_pos] == "+" or self.bases[base_pos] == "-":
                #Remove insertion
                indelLength = self.__getIndelLength(self.bases,base_pos)
                base_pos += (indelLength + int(len(str(indelLength))) +1)
                continue
            
            if self.bases[base_pos] == "^":
                #Remove '^x' character
                base_pos += 2
                continue

            if self.bases[base_pos] in validBases:
                #valid base
                self.bases_masked = self.bases_masked + self.bases[base_pos]
                base_pos += 1
                continue

            else:
                base_pos +=1

        if len(self.baseQualities) != len(self.bases_masked):
            print "basemask"
            #print str(self.bases)
            #print str(self.bases_masked)
            #print str(self.baseQualities)
        return
    
    def __getBasePositions(self):
        """
        Counts the position in self.bases_masked of each base type.
        Storing the info in self.basePositions
        """
        for i in range(len(self.bases_masked)):
            char = self.bases_masked[i]        
            self.basePositions[char] = self.basePositions.get(char,[])
            self.basePositions[char].append(i)
        
        self.refBaseCount = len(self.basePositions["."] + self.basePositions[","])
        
        return 

    def __getMajorAltBase(self):
        """
        Determines the major alt base
        """
        countDict = {"Aa":len(self.basePositions["A"] + self.basePositions["a"]),\
            "Tt":len(self.basePositions["T"] + self.basePositions["t"]),\
            "Cc":len(self.basePositions["C"] + self.basePositions["c"]),\
            "Gg":len(self.basePositions["G"] + self.basePositions["g"])}
       
        maxItem = max(countDict.iteritems(),key=operator.itemgetter(1))
        if int(maxItem[1]) > 0:
            self.majorAltBase = maxItem[0]
            self.majorAltBaseCount = int(maxItem[1])
        
        return 

    def __repr__(self):
        """
        Return a string representation of this mpSample object
        """
        repr_str = "bases: " + str(self.bases) + "\n" + \
             "bases_masked: " + str(self.bases_masked) + "\n" + \
             "baseQualities : " + str(self.baseQualities) + "\n" + \
             "depth: " + str(self.depth) + "\n" + \
             "mappingQual: " + str(self.mappingQual) + "\n" + \
             "majorAltBase: " + str(self.majorAltBase) + "\n" + \
             "majorAltBaseCount: " + str(self.majorAltBaseCount) + "\n"

        return repr_str

    def altReadCount(self,min_qual=0):# --> int
        """
        Returns the number of major alternate reads present with base quality > min_qual
        """
        good_base = 0
        for base_pos in self.basePositions[self.majorAltBase[0]] + self.basePositions[self.majorAltBase[1]]:
            if (ord(self.baseQualities[base_pos]) - 33) >= min_qual:
                good_base += 1

        return good_base



    def altReadBothStrands(self,min_qual=0):# -- bool
        """
        Returns whether there is a major alt read on both strands
        """
        forward_base = self.basePositions[self.majorAltBase[0]] 
        reverse_base = self.basePositions[self.majorAltBase[1]]

        if len(forward_base) > 0 and len(reverse_base) > 0:
            #Atleast one in each direction
            if any((ord(self.baseQualities[base_index]) - 33) >= min_qual for base_index in forward_base):
                if any((ord(self.baseQualities[base_index]) - 33) >= min_qual for base_index in reverse_base):
                    return True

        return False


    def getBaseCount(self,base):# --> int
        """
        Return the number of reads for a specific base.
        DIRECTION OF BASES MATTERS A vs. a
        """
        return len(self.basePositions.get(base,[]))




