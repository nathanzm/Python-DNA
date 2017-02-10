import BioDNA # allows for use of defined functions within BioDNA.py

def main():
    
    print ("\n         ---+++ GeneAnalysis (build v3.02.14) +++---")
    print ("                     - by Nathan Morse\n")
    AAtable = BioDNA.makeAminoAcidTable() # instance of amino acid table
    
    #for (tripleRNA, AAcode) in AAtable.items():
    #    print (tripleRNA, AAcode)

        
    DNA = "ATGTGTACGCAAATGATATCGTATTAG"; # string of DNA for input
    
    antiSenseStrandDNA = BioDNA.complementDNA_to_DNA(DNA) # complementary DNA
    mRNA = BioDNA.transcribe(DNA) # stores call to function transcribe on DNA
    RNAlength = len(mRNA) # the length of the RNA

    print (" DNA: 5'", DNA, "3'") # prints the DNA strand
    print (" DNA: 3'", antiSenseStrandDNA, "5'") # prints the anti-sense strand
    
    # --------- TRANSCRIBE --------------------------

    print ("        ", len(DNA) * "|") # prints according number of match lines
    print ("mRNA: 5'", mRNA, "3'") # prints the mRNA
      

    print ("\n==============================================================")
    print ("Reading Frame #1")
    # --------- TRANSLATE Reading Frame #1 --------------------------
    
    proteins = BioDNA.translate(mRNA, AAtable) # translates mRNA to proteins
    numberOfconnectors = int((len(proteins)) / 3) # generates num of match lines

    print ("RNA length:", RNAlength, "bp") # prints length of mRNA
    print ("      RNA:        ", mRNA) # prints the mRNA
    print ("                  ", numberOfconnectors * " | ") # prints connectors
    print ("      Protein:    ", proteins, "\n") # prints protein chain
    
    BioDNA.readingFrameChecks(proteins) # prints information about amino acids
    
    print ("==============================================================\n")
    
    
    print ("\n==============================================================")
    print ("Reading Frame #2")
    # --------- TRANSLATE Reading Frame #2 --------------------------
    """
    reading frame #2 does exactly as the previous reading frame except it starts
    the analysis of the mRNA 1 position further
    """
    proteins = BioDNA.translate(mRNA[1:], AAtable)
    numberOfconnectors = int((len(proteins)) / 3)

    print ("RNA length:", RNAlength - 1, "bp") # subtract 1 from length
    print ("      RNA:        ", mRNA[1:])
    print ("                  ", numberOfconnectors * " | ")
    print ("      Protein:    ", proteins, "\n")
    
    BioDNA.readingFrameChecks(proteins)
    
    print ("==============================================================\n")
    
    
    print ("\n==============================================================")
    print ("Reading Frame #3")
    # --------- TRANSLATE Reading Frame #3 --------------------------
    """
    reading frame #3 does exactly as the previous reading frame except it starts
    the analysis of the mRNA 1 position further
    """
 
    proteins = BioDNA.translate(mRNA[2:], AAtable)
    numberOfconnectors = int((len(proteins)) / 3)

    print ("RNA length:", RNAlength - 2, "bp") # subtract 2 from length
    print ("      RNA:        ", mRNA[2:])
    print ("                  ", numberOfconnectors * " | ")
    print ("      Protein:    ", proteins, "\n")
    
    BioDNA.readingFrameChecks(proteins) 
    
    print ("==============================================================\n")

# --- end main() ---------


#---------------------------------------------------------
# Python starts here ("call" the main() function at start
if __name__ == '__main__':
    main()
#---------------------------------------------------------