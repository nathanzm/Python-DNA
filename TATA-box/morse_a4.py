"""----------------------------------------------------------------------------
SUMMARY: This program investigates the DNA sequences upstream to a given gene
based on the location of what is known as the TATA-box. The TATA-box consists
of the nucleotides TATA followed by either an A or T, then an A, then an A or T,
then either an A or G. The program finds various motifs through use of regular
expressions within the upstream region of the gene prior to the TATA-box (provided
that both exist, but the program will test for existence and prevents any error
if there is no TATA-box or upstream/downstream region).

Programmer: Nathan Morse
    
INPUT: A FASTA formatted input file with the DNA of an organism; could be
a specific gene.

OUTPUT: The size of the entire DNA strand given in base pairs along with the actual
DNA given lined up with its according index numbers. The location of the upstream, 
downstream, and TATA-box are shown accordingly as well as the actual upstream of 
the TATA-box (aligned with index numbers) and its length in base pairs. Any direct
repeat motifs of the upstream region are printed with their locations as well as
percent of direct repeats in the upstream. Any mirror repeat motifs in the upstream
region are printed with their locations and the percentage of mirror motifs of
the upstream region is printed also. The favorite motif of A or T followed by 
another A or T is found within the upstream and printed with its according location
and percent of the upstream region as well. If any of the above mentioned do not
exist, an error message will be printed appropriately.
-------------------------------------------------------------------------------
"""

import re, os # import regex library and os for checking file inputs
# ------------------------------------------------------------------------

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
#--------------------(end of printIndexing(string)-------------------------

def TATAbox(retrievedDNA):
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
            
            # call to regex functions on upstream:
            if ("upstream" in locals()): # if upstream variable exists locally
                findDirectRepeats(upstream)
                findMirrorRepeats(upstream)
                findATRepeats(upstream)
        else: # outtermost else statement in case there is no TATA-box:
            print("There is no TATA box to locate.")
#-------------------(end of TATAbox())-------------------------------------

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
#-----------------(end of findDirectRepeats(upstream))---------------------

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
#-----------------(end of findMirrorRepeats(upstream))---------------------

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
#-----------------(end of findCATRepeats(upstream))---------------------

def main():
    
    inputFile = input("Enter a file name: ") # user input file name
    retrievedDNA = getDNA(inputFile) # call to getDNA on inputted file
    TATAbox(retrievedDNA) # call to TATAbox function on the retrievedDNA.
    

#---------------------(end of main)-----------------------
# Python starts here ("call" the main() function at start)
if __name__ == '__main__':
    main()
#---------------------------------------------------------