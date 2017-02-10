"""BioDNA : a module of DNA processing functions """

 # libraries for use of regex, glob, counter, and checking file inputs:
import re, glob, os
from collections import Counter

# ******* THESE FUNCTIONS WERE ORIGINALLY MADE FOR PROGRAM ASSIGNMENT 3 ********

#--------------------------------------------------------------------------
def complementDNA_to_DNA( DNA ):
    """This function returns the complement strand of DNA.
    -----------------------------------------------------------------------    
    """

    directStrandTransfer = "ACTG" # transform these following
    inDirectStrandTransfer = "TGAC" # into these accordingly:
    transTable = str.maketrans(directStrandTransfer, inDirectStrandTransfer)
    complementaryStrand = DNA.translate(transTable) # stores anti-sense strand
    
    return complementaryStrand
#----------------------(end of complementDNA_to_DNA())---------------------

#--------------------------------------------------------------------------
def complementDNA_to_RNA( DNA ):
    """ This function returns the complement strand of RNA using the same
    method as above.
    -----------------------------------------------------------------------
    """

    directStrandTransfer = "ACTG"
    RNAinDirectStrandTransfer = "UGAC"
    transTable = str.maketrans(directStrandTransfer, RNAinDirectStrandTransfer)
    complementaryStrand = DNA.translate(transTable)

    return complementaryStrand
# ----------------------(end of complementDNA_to_RNA())--------------------

#--------------------------------------------------------------------------
def transcribe( DNA ):
    """ This function simulates transcription by building template strand and
    then complementing that to RNA.
    -----------------------------------------------------------------------
    """
    # gene is on the DNA strand
    # but RNA copy is built off the template (or anti-sense) strand

    mRNA = "" # create empty string to put mRNA in

    # (a) get complement of the DNA strand (that is, make the template strand)
    
    DNAtemplateStrand = complementDNA_to_DNA(DNA)
    
    # (b) get mRNA complement of the template strand
    
    mRNA += complementDNA_to_RNA(DNAtemplateStrand) # put into empty string (no longer empty)
    
    # (c) return the mRNA
    
    return mRNA
#------------------------(end of transcribe())-----------------------------

#--------------------------------------------------------------------------
def makeAminoAcidTable( ):
    """ Returns a Python Dictionary (hash table) that maps RNA nucleotides to Amino Acid symbols;
    both one letter and three letter symbols are returned; the user will need to parse which of
    these symbols they want to use.
    -----------------------------------------------------------------------    
    """
    AAtable = { "UUU":"F|Phe","UUC":"F|Phe","UUA":"L|Leu","UUG":"L|Leu","UCU":"S|Ser","UCC":"S|Ser",
                "UCA":"S|Ser","UCG":"S|Ser","UAU":"Y|Tyr","UAC":"Y|Tyr","UAA":"*|***","UAG":"*|***",
                "UGU":"C|Cys","UGC":"C|Cys","UGA":"*|***","UGG":"W|Trp","CUU":"L|Leu","CUC":"L|Leu",
                "CUA":"L|Leu","CUG":"L|Leu","CCU":"P|Pro","CCC":"P|Pro","CCA":"P|Pro","CCG":"P|Pro",
                "CAU":"H|His","CAC":"H|His","CAA":"Q|Gln","CAG":"Q|Gln","CGU":"R|Arg","CGC":"R|Arg",
                "CGA":"R|Arg","CGG":"R|Arg","AUU":"I|Ile","AUC":"I|Ile","AUA":"I|Ile","AUG":"M|Met",
                "ACU":"T|Thr","ACC":"T|Thr","ACA":"T|Thr","ACG":"T|Thr","AAU":"N|Asn","AAC":"N|Asn",
                "AAA":"K|Lys","AAG":"K|Lys","AGU":"S|Ser","AGC":"S|Ser","AGA":"R|Arg","AGG":"R|Arg",
                "GUU":"V|Val","GUC":"V|Val","GUA":"V|Val","GUG":"V|Val","GCU":"A|Ala","GCC":"A|Ala",
                "GCA":"A|Ala","GCG":"A|Ala","GAU":"D|Asp","GAC":"D|Asp","GAA":"E|Glu",
                "GAG":"E|Glu","GGU":"G|Gly","GGC":"G|Gly","GGA":"G|Gly","GGG":"G|Gly"}
    
    return AAtable
# ------------------------(end of makeAminoAcidTable())-------------------- 

#--------------------------------------------------------------------------
def translate( mRNA, AAtable ):  
    """ This function returns the 3-letter amino acid string of symbols for 
    the corresponding mRNA.
    -----------------------------------------------------------------------
    """
    protein = ""
    
    startBP = 0
    endBP   = len(mRNA)
    
    nextCodonStart = startBP
    # while we don't run off the end of the sequence of mRNA, snag nucleotide codons
    #       and build a sequence of amino acid symbols (uses the 3-letter symbols)
    while ( nextCodonStart+3 <= endBP):
        
        # slice out the new RNA nucleotide triple
        nextCodon = mRNA[nextCodonStart:nextCodonStart+3]
        
        # play RIBOSOME: convert an RNA nucleotide-triple to an amino acid symbol
        nextAA = AAtable[ nextCodon ]
        
        # amino acid symbols come as two types of symbols (3 letter and 1 letter);
        # in general:  triple:<one letter code>|<three letter code>
        # e.g.,  "GAG":"E|Glu" ... thus we need to parse out which type of symbol we want
        
        # we want the 3-letter 'Glu' version here, so slice out from location 2 to the end
        AAsymbol = nextAA[2:]
        
        # add this AA onto the end of the growing amino-acid (protein) chain
        protein = protein + AAsymbol
        
        # move down to next triple
        nextCodonStart = nextCodonStart + 3
        
    # end while more triples to check
    
    return protein
# -------------------------(end of translate())----------------------------

#--------------------------------------------------------------------------
def readingFrameChecks ( proteins ):
    """ This function simply checks conditions of start/stop amino acids and
    generates output for reading frames
    -----------------------------------------------------------------------
    """
    if (proteins[:3] == "Met"): # if there is a Met at the first 3 characters
        print ("There is a Met at the start site.")
    else:
        print ("There is no Met at the start.")
    if (proteins[-3:] == "***"): # if there is stop amino acid at the last 3 characters
        print ("There is a stop at the end.")
    else:
        print ("There is no stop at the end.")
    if (proteins.find("Met", 3) != -1): # if there is Met somehwere else after the start
        print ("There is at least one Met somewhere other than start.")
    else:
        print ("There are no Mets in positions outside of the start site.")
    # and lastly if there is not a stop amino acid at the end and there is one elsewhere:
    if (proteins.find("***") != (len(proteins) - 3) and proteins.find("***") != -1):
        print ("There is at least one STOP interrupting the gene.")
    else:
        print ("There are no interrupting STOPS.")
#-------------------(end of readingFrameChecks())--------------------------

# ********* END OF FUNCTIONS ORIGINALLY MADE FOR PROGRAM ASSIGNMENT 3 **********

# ******* THESE FUNCTIONS WERE ORIGINALLY MADE FOR PROGRAM ASSIGNMENT 4 ********

#--------------------------------------------------------------------------
def getDNA(filename):
    """ Open a FASTA file of DNA read it, and return the DNA as one string.
    
    Function to open a FASTA formatted file of DNA sequence, remove
    first (header) line and newline characters, and return as one (long) string.
    
    Argument:  one(1) string with the name of a FASTA-formtted file of DNA.
    Returns:   entire sequence of DNA as one string
 
    A sample of how you might "call" this funciton from your program
    to 'get' the DNA for chromosome III of the worm, C.elegans.
 
    sequence = getDNA("worm_III.fna")

    Note: This function (sometimes called a subroutine) is ready for you
    to use in your program. Think of this routine as a "button on your
    calculator" ... this is ready for to use. When you "use" this
    routine, we say "your program is CALLING" for the use of that routine.
    -----------------------------------------------------------------------
    """
    # make sure there really *is* a filename with this name
    if (not os.path.isfile(filename)):
        print("No file found in current directory named: ", filename)
        return ""
    else:
        DNA = ""
        with open(filename) as INPUT:
            next(INPUT)  # skip the header line
            # now read all the other lines in the file
            for nextLine in INPUT: 
                nextLine = nextLine.strip()  # remove the newline
                DNA = DNA + nextLine
            # end for each line
                
            return DNA.lower() # make DNA lowercase
    # end else
#-------------------(end of getDNA())--------------------------------------

#--------------------------------------------------------------------------
def TATAboxAnalysis(retrievedDNA):
    """ This function finds the first TATA box it finds starting from the front
    of the sequence as well as executes all other functions of the program accordingly 
    (other than getDNA(filename)) based on if there is an upstream region to the
    TATA-box or not.
    -----------------------------------------------------------------------
    """   
    if (retrievedDNA != ""): # if there is DNA   
        TATAregex = re.compile(r"tata[at]a[at][ag]") # regex for TATA-box
    
        if (TATAregex.search(retrievedDNA)): # if the TATA-box exists:
            TATAstart = (TATAregex.search(retrievedDNA)).start() # start of TATA-box
            TATAend = (TATAregex.search(retrievedDNA)).end() # end of TATA-box
            # the next line puts only the TATA-box in uppercase (all in a stored variable):
            retrievedDNA = retrievedDNA[:TATAstart] + retrievedDNA[TATAstart:TATAend].upper() + retrievedDNA[TATAend:]
            
            print("\nThe size of the region being searched is", len(retrievedDNA), "bp")
            print("=========================================================")
            print(retrievedDNA) # print the sequence with the TATA-box in uppercase
            
            printIndexing(retrievedDNA) # call to printIndexing function on retrievedDNA
                    
            if (TATAstart != 0): # if there is an upstream region
                print("\n\nUpstream:   [ 1 :", TATAstart, "]") # print location
                upstream = retrievedDNA[:TATAstart] # store location of upstream
                ''' 
                the next 2 lines calculate the percentage of the upstream
                region and make that public for use in outside functions
                
                '''
                global percentageUpstream
                percentageUpstream = (len(upstream)/len(retrievedDNA)*100)
            else: # otherwise there is no upstream region to the TATA-box!
                print("\n\nThere is no upstream region prior to TATA-box!")
            
            # now print the location (start to end) of TATA-box    
            print("TATA-box:   [", TATAstart+1, ":", TATAend, "]")
            
            if (TATAend != len(retrievedDNA)): # if there is a downstream, print
                print("Downstream: [", TATAend+1, ":", len(retrievedDNA), "]")
            else: # if the end of TATA-box is at the end of the sequence also:
                print("There is no downstream region after the TATA-box!")
            
            print("---------------------------------------------------------")
            
            # call to regex functions on upstream (these are defined below):
            if ("upstream" in locals()): # if upstream variable exists locally
                findDirectRepeats(upstream)
                findMirrorRepeats(upstream)
                findATRepeats(upstream)
        else: # outtermost else statement in case there is no TATA-box:
            print("There is no TATA box to locate.")
#-------------------(end of TATAboxAnalysis())-----------------------------

#--------------------------------------------------------------------------
def printIndexing(string):
    """ This function prints the indexing (by tens) based on the length of the
    region to help the reader verify locations (bp) of motifs.
    -----------------------------------------------------------------------
    """    
    nucleotideCount = 1 # start reading position 1 rather than 0
    
    for nucleotide in string: # for each A, T, G, or C
        if (nucleotideCount == 10): # if count gets to 10,
            nucleotideCount = 0 # set back to 0
        if (nucleotideCount != 10):
            print(nucleotideCount, end="") # don't print on a new line each time
            nucleotideCount += 1 # advance the count
    # end of for loop
#--------------------(end of printIndexing()-------------------------------

#--------------------------------------------------------------------------
def findDirectRepeats(upstream):
    """ This function finds all direct repeat motifs with a total length 
    between 4 to 12 bp and their location within the substring. It also reports
    the percentage of the sum of all direct repeats as compared to the entire 
    substring (to the left of the TATA-box).
    -----------------------------------------------------------------------
    """       
    print("\n\nUPSTREAM of TATA-box")
    print("---------------------------------------------------------")
    print(upstream) # print everything before the TATA-box
    printIndexing(upstream) # call to printIndexing on upstream
    print("\n\nThe size of the upstream region is", len(upstream), "bp")
    # the next line prints the percentage of the upstream region by accessing the public variable
    print("The percentage of upstream region is:{:5.1f}%".format(percentageUpstream))
    print("---------------------------------------------------------")
    print("-Searching for Direct Repeats (DR) in the upstream region-\n")
    
    DRregex = re.compile(r"(.{2,6})\1") # regex for finding direct repeats
    DRiterator = DRregex.finditer(upstream) # find all matches
    
    DRlengths = 0 # counter for lengths of direct repeats
    
    for nextDR in DRiterator: # loop through list of matches
        DR = nextDR.group() # get the actual regex match
        DRstart = nextDR.start() # where match begins
        DRend = nextDR.end() # where match ends
        DRlengths = DRlengths + (DRend - DRstart) # length of current match added to DRlengths
        print("Found DR:", DR)
        print("     at upstream location: [", DRstart+1, ":", DRend, "]")
    # end of for loop
        
    percentageDR = (DRlengths / len(upstream)) * 100 # calculate % of direct repeat
    print("\n\nPercent of DR's in the upstream region is:{:5.1f}%".format(percentageDR))
    print("---------------------------------------------------------")
    print("---------------------------------------------------------")
#-----------------(end of findDirectRepeats())-----------------------------

#--------------------------------------------------------------------------
def findMirrorRepeats(upstream):
    """ This function works exactly the same way as the previous 
    (findDirectRepeats(upstream)) except it reports all mirror repeats of
    length 6 bp and their location within the substring.
    -----------------------------------------------------------------------
    """ 
    print("-Searching for Mirror Repeats (MR) in the upstream region-\n")
        
    MRregex = re.compile(r"(.)(.)(.)\3\2\1") # regex for mirror repeats of length 6
    MRiterator = MRregex.finditer(upstream)
        
    MRlengths = 0
        
    for nextMR in MRiterator:
        MR = nextMR.group()
        MRstart = nextMR.start()
        MRend = nextMR.end()
        MRlengths = MRlengths + (MRend - MRstart)
        print("Found MR:", MR)
        print("     at upstream location: [", MRstart+1, ":", MRend, "]")
    # end of for loop
            
    percentageMR = (MRlengths / len(upstream)) * 100
    print("\n\nPercent of MR's in the upstream region is:{:5.1f}%".format(percentageMR))
    print("---------------------------------------------------------")
    print("---------------------------------------------------------")    
#-----------------(end of findMirrorRepeats())-----------------------------

#--------------------------------------------------------------------------
def findATRepeats(upstream):
    """ This function works exactly the same way as the previous two except it 
    reports all occurences of a favorite motif and their location within the
    substring. The motif is any a or t followed by another a or t, repeated
    1 to 8 times.
    -----------------------------------------------------------------------
    """ 
    print("-Searching for AT/TA/AA/TT runs in the upstream region-\n")
        
    ATregex = re.compile(r"([at][at]){1,8}") # regex for favorite motif as described
    ATiterator = ATregex.finditer(upstream)
        
    ATlengths = 0
        
    for nextAT in ATiterator:
        AT = nextAT.group()
        ATstart = nextAT.start()
        ATend = nextAT.end()
        ATlengths = ATlengths + (ATend - ATstart)
        print("Found favorite motif:", AT)
        print("     at upstream location: [", ATstart+1, ":", ATend, "]")
    # end of for loop
            
    percentageAT = (ATlengths / len(upstream)) * 100
    print("\n\nPercent of AT/TA/AA/TT runs in the upstream region is:{:5.1f}%".format(percentageAT))
    print("---------------------------------------------------------")
#-----------------(end of findCATRepeats())--------------------------------

# ********* END OF FUNCTIONS ORIGINALLY MADE FOR PROGRAM ASSIGNMENT 4 **********

# ******* THESE FUNCTIONS WERE ORIGINALLY MADE FOR PROGRAM ASSIGNMENT 5 ********

#--------------------------------------------------------------------------
def inputDirectoryLengths(inputDirectory):
    """ The following function goes through a given input directory and prints
    the length of DNA found in each file (using the getDNA()) function in 
    addition to printing the name of each file. All non-FASTA formatted files
    are ignored.
    -----------------------------------------------------------------------
    """
    if (os.path.isdir(inputDirectory)): # if the directory does exist
        fullPath = os.path.join(inputDirectory, "*.fna") # only get .fna files

        # put all FASTA files (*.fna) in the test directory in a list called fileList
        fileList = glob.glob(fullPath)
        
        if (fileList == []): # if the fileList is empty
            print("The directory does not contain any FASTA formatted files!")
        else:
            fileList.sort() # put the file names in alphabetical order; easier to see
            print("The files within", inputDirectory, "are as follows:")

            for nextFile in fileList:
                if (os.path.isfile(nextFile)):
                    # ok, it is SAFE to open the file
                    INPUT = open(nextFile, 'r')
                    currentFile_DNALength = len(getDNA(nextFile)) # length of DNA sequence in file
                    print(nextFile, "with DNA length:", currentFile_DNALength, "bp.")
                    
                    INPUT.close() # close input file
                    
                else:
                    print("Sorry, there is no file called: ", nextFile)
            # end for each file in fileList
    else: # the directory does not exist
        print("The directory", inputDirectory, "does not exist!")
#-----------------------(end of inputDirectoryLengths())-------------------

#--------------------------------------------------------------------------
def breakIntoMotifs(inputDirectory, LmerSize, numberOfTopRankings):
    """ The following function goes through a given input directory and its files
    and uses the LmerSize given through user input as well as the numberOfTopRankings
    to write the file name, total number of motifs, the given number of top ranking
    motifs with their motif name, frequency, and propotion to an output file (.csv).
    -----------------------------------------------------------------------
    """    
    motifDictionary = {} # dictionary of motifs with number of occurrences
    AllLmers = [] # list containing ALL motifs based on LmerSize given
    
    if (os.path.isdir(inputDirectory)): # if the directory does exist
        fullPath = os.path.join(inputDirectory, "*.fna") # only get .fna files

        # put all FASTA files (*.fna) in the test directory in a list called fileList
        fileList = glob.glob(fullPath)
        
        if (fileList != []): # if the fileList is not empty
            fileList.sort() # put the file names in alphabetical order; easier to see

            OUTPUT = open("Results.csv", 'w') # open file for output
            COMMA = "," # define COMMMA for comma separated values

            for nextFile in fileList: # for each file in the fileList
                if (os.path.isfile(nextFile)):
                    # ok, it is SAFE to open the file
                    INPUT = open(nextFile, 'r')
                    DNASequence = getDNA(nextFile) # use getDNA to put DNA sequence from file in variable
                    startOfmotif = 0 # start of motif starts at the first nucleotide in DNA strand
                    endOfmotif = LmerSize # end of motif is the size of the L-mer
                    Lmer = DNASequence[startOfmotif:endOfmotif] # the first motif in DNA strand
                    numberOfmotifs = len(DNASequence) - (LmerSize - 1) # the total number of motifs
                    
                    INPUT.close() # close input file; no longer needs to be open
                    
                    iterator = 0 # used for keeping track of shifting the motif in the DNA strand
                    while (iterator != numberOfmotifs): # while the iterator isn't at the last motif
                        AllLmers.append(Lmer.upper()) # add the uppercase motif to AllLmers list
                        iterator += 1 # advance the iterator
                        startOfmotif += 1 # shift the start position of the motif over one
                        endOfmotif += 1 # shift the end position of the motif over one
                        Lmer = DNASequence[startOfmotif:endOfmotif] # the L-mer is shifted
                    
                    for motif in AllLmers: # for each motif in AllLmers
                        if (motif in motifDictionary): # if the motif is already in the dictionary:
                            motifDictionary[motif] = motifDictionary[motif] + 1 # increment that motif's count
                        else: # first time this motif is spotted
                            motifDictionary[motif] = 1 # the count for that motif is 1
                    # end for each motif in AllLmers
                    
                    """ The next line "simulates" an ordered dictionary by using the built in function Counter from
                    the collections library to collect the counts from the values in the motifDictionary and then
                    uses the built in function .most_common to get the given top number of counts (indicated by the
                    variable numberOfTopRankings) from the collected counts. By using this method, the top number
                    of motif occurrences (values in the motifDictionary) along with their according motif names 
                    (keys in the motifDictionary) can be put into the variable orderedDictionary. The variable
                    orderedDictionary is not actually a dictionary but is of type list.
                    """
                    orderedDictionary = Counter(motifDictionary).most_common(numberOfTopRankings)
                        
                    # write headers for output in output file (all cells are separated using commas)
                    OUTPUT.write("%s%s\n%s%i\n" % ("FILE: ", nextFile, "Number of motifs: ", numberOfmotifs))
                    OUTPUT.write("%s%c%s%c%s%c%s\n" % ("Rank:", COMMA, "Motif:", COMMA, "Frequency:", COMMA, "Proportion:"))
                    
                    # write dictionary contents in output file
                    currentRanking = 1 # counter used to keep track of printing ranking
                    for k, v in orderedDictionary: # for the "keys" and "values" in orderedDictionary
                        motifProportion = v / numberOfmotifs # calculate proportion of motif
                        OUTPUT.write("%i%c%s%c%i%c%f\n" % (currentRanking, COMMA, k, COMMA, v, COMMA, motifProportion))
                        currentRanking += 1 # advance the current ranking
                            
                    OUTPUT.write("\n") # separate each file with a newline
                    
                    motifDictionary.clear() # clear the motif dictionary for use in the next file
                    del AllLmers[:] # clear the AllLmers list for use in the next file
                    
                # end of for nextFile in fileList                       
                    
                else:
                    print("Sorry, there is no file called: ", nextFile)
            # end for each file in fileList
        
        OUTPUT.close() # close the output file 
        print("The output file has been successfully created!")         
        
    else: # the file directory given does not exist
        return  
#-----------------------(end of breakIntoMotifs())------------------------

# ********* END OF FUNCTIONS ORIGINALLY MADE FOR PROGRAM ASSIGNMENT 5 **********