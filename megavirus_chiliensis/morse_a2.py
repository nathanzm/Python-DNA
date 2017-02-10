import os
from math import exp, log

"""----------------------------------------------------------------------------
SUMMARY: This python program reads a FASTA formatted input file filled with DNA
from a given organism and prints a neat summary of "Chargaff's numbers", which 
contains a helpful analysis of the DNA.
    
Programmer: Nathan Morse
    
INPUT: A FASTA formatted input file filled with the DNA of an organism.

OUTPUT: The total number of each of the different nucleotides along with their 
percentages and proportions within the DNA, the total number of nucleotides in 
the entire file, the total number of nucleotides in the file that aren't A, C, 
G, or T. The output also includes the likelihood of finding the specific 8-mer
motif ACGTACGT as well as the length in centimeters and meters of a strand of 
DNA equal to 100,000 bp.
-------------------------------------------------------------------------------
"""

 #------- Do NOT change this Python function--------------------------------

def getDNA( filename ):
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
        print ("No file found in current directory named: ", filename)
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
                
            return DNA
    # end else
    
# end getDNA()
# ------------------------------------------------------------------------
    
    
def main():    
    
    inputGenomeFile = "megavirus_chiliensis.fna" # file for input
    # or
    # inputGenomeFile = "Nanoarchaeum_equitans_Kin4_M.fna"
    
    # get the DNA sequence from a file and store in a Python string
    sequence = getDNA(inputGenomeFile)
     
    # don't run this next line if you have a huge genome!
    # print ("\nDNA: ", sequence, "\n") # prints the DNA sequence
    
    DNA = sequence.upper() # makes sure the sequence is capitalized
    DNAlength = len(DNA) # the length of the DNA sequence
    print("Total length of DNA:\t", DNAlength, "bp") # prints DNAlength
    print("-------------------------------------------------------------------------")
    #-----------------------------------------
    
    # COUNT NUCLEOTIDES:    
    # first store the count of each nucleotide in a variable
    adenineCount = DNA.count("A") # number of adenine nucleotides
    cytosineCount = DNA.count("C") # number of cytosine nucleotides
    guanineCount = DNA.count("G") # number of guanine nucleotides
    thymineCount = DNA.count("T") # number of thymine nucleotides
    
    # total number of nucleotides that aren't A, C, G, or T
    otherCount = 0 # used to keep track of nucleotides that are not A, C, G, or T
    for i in DNA: # for all the characters in DNA:
        if (i != "A" and i != "C" and i != "G" and i != "T"):
            otherCount += 1 # advance otherCount by 1 if there is something other than A, C, G, or T
    
    # now print out the formatted number of each nucleotide (stored in the dedicated variables)
    print("The number of Adenine (A) nucleotides in the file of DNA:\t{:5}".format(adenineCount))
    print("The number of Cytosine (C) nucleotides in the file of DNA:\t{:5}".format(cytosineCount))
    print("The number of Guanine (G) nucleotides in the file of DNA:\t{:5}".format(guanineCount))
    print("The number of Thymine (T) nucleotides in the file of DNA:\t{:5}".format(thymineCount))
    
    # now print out the number of nucleotides that aren't A, C, G, or T
    print("The number of nucleotides that are not A, C, G, or T:\t{:5}".format(otherCount))
    print("-------------------------------------------------------------------------")
    #-----------------------------------------
    
    # PROPORTIONS OF NUCLEOTIDES:
    # first store the proportions in variables; these are calculated by count divided by total length
    adenineProportion = adenineCount/DNAlength
    cytosineProportion = cytosineCount/DNAlength
    guanineProportion = guanineCount/DNAlength
    thymineProportion = thymineCount/DNAlength
    
    # now print the formatted proportions of each nucleotide
    print("The proportion of Adenine (A) nucleotides in the file of DNA:\t{:5.2f}".format(adenineProportion))
    print("The proportion of Cytosine (C) nucleotides in the file of DNA:\t{:5.2f}".format(cytosineProportion))
    print("The proportion of Guanine (G) nucleotides in the file of DNA:\t{:5.2f}".format(guanineProportion))
    print("The proportion of Thymine (T) nucleotides in the file of DNA:\t{:5.2f}".format(thymineProportion))
    print("-------------------------------------------------------------------------")
    #-----------------------------------------
    
    # PERCENTAGES OF PROPORTIONS:
    # first store the percentages in variables; these are calculated by proportion times 100
    adeninePercentage = adenineProportion * 100
    cytosinePercentage = cytosineProportion * 100
    guaninePercentage = guanineProportion * 100
    thyminePercentage = thymineProportion * 100
    
    # now print the formatted percentages of each nucleotide
    print("The percentage of Adenine (A) nucleotides in the file of DNA:\t{:5.2f}".format(adeninePercentage), "%")
    print("The percentage of Cytosine (C) nucleotides in the file of DNA:\t{:5.2f}".format(cytosinePercentage), "%")
    print("The percentage of Guanine (G) nucleotides in the file of DNA:\t{:5.2f}".format(guaninePercentage), "%")
    print("The percentage of Thymine (T) nucleotides in the file of DNA:\t{:5.2f}".format(thyminePercentage), "%")
    print("-------------------------------------------------------------------------")
    #-----------------------------------------------------------
    
    # compute liklihood of ACGTACGT (8-mer) in this sequence:
    motifLog = log(adenineProportion * cytosineProportion * guanineProportion * thymineProportion * adenineProportion * cytosineProportion * guanineProportion * thymineProportion)
    # take the log of all proportions of ACGTACGT multiplied together then calculate:
    motifProbability = exp(motifLog) # e^motifLog
    
    print("The liklihood of ACGTACGT (8-mer) in this sequence:\t{:5.2f}".format(motifProbability))
    print("-------------------------------------------------------------------------")
    #-----------------------------------------------------------------------

    # compute the length of a piece of DNA in centimeters (cm):
    sampleLengthBP = 100000 # length of strand of DNA
    mmLength = sampleLengthBP/4000000 # there are 4,000,000 bases in a millimeter of DNA
    cmLength = mmLength/10 # converts millimeters to centimeters
    mLength = cmLength/100 # converts centimeters to meters
    
    # results printed in scientific notation
    print("The length (in cm) of a strand of DNA equal to 100,000 bp:\t{:5.2e}".format(cmLength))
    print("The length (in m) of a strand of DNA equal to 100,000 bp:\t{:5.2e}".format(mLength))
    print("-------------------------------------------------------------------------")

# === end main() ==========================================
    
    
    
#-----------------------------------------------------
# Python starts here
if __name__ == '__main__':
    main()
#-----------------------------------------------------  