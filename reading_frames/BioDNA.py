"""BioDNA : a module of DNA processing functions """

# ----------------------------------------
def complementDNA_to_DNA( DNA ):
    """Returns the complement strand of DNA"""

    directStrandTransfer = "ACTG" # transform these following
    inDirectStrandTransfer = "TGAC" # into these accordingly:
    transTable = str.maketrans(directStrandTransfer, inDirectStrandTransfer)
    complementaryStrand = DNA.translate(transTable) # stores anti-sense strand
    
    return complementaryStrand
# ----------------------------------------

# ----------------------------------------
def complementDNA_to_RNA( DNA ):
    """Returns the complement strand of RNA using the same method as above"""

    directStrandTransfer = "ACTG"
    RNAinDirectStrandTransfer = "UGAC"
    transTable = str.maketrans(directStrandTransfer, RNAinDirectStrandTransfer)
    complementaryStrand = DNA.translate(transTable)

    return complementaryStrand
# ----------------------------------------

# ----------------------------------------    
def transcribe( DNA ):
    """Simulates transcription by building template strand and then complementing that to RNA"""
    
    # gene is on the DNA strand
    # but RNA copy is built off the template (or anti-sense) strand

    mRNA = "" # create empty string to put mRNA in

    # (a) get complement of the DNA strand (that is, make the template strand)
    
    DNAtemplateStrand = complementDNA_to_DNA(DNA)
    
    # (b) get mRNA complement of the template strand
    
    mRNA += complementDNA_to_RNA(DNAtemplateStrand) # put into empty string (no longer empty)
    
    # (c) return the mRNA
    
    return mRNA
# ---------------------------------------- 

# ---------------------------------------- 
    
def makeAminoAcidTable( ):
    """Returns a Python Dictionary (hash table) that maps RNA nucleotides to Amino Acid symbols;
    both one letter and three letter symbols are returned; the user will need to parse which of
    these symbols they want to use"""
    
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
# ---------------------------------------- 

# ---------------------------------------- 

def translate( mRNA, AAtable ):  
    """Returns the 3-letter amino acid string of symbols for the corresponding mRNA"""
    
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
# ---- end translate() -------------------

# ---------------------------------------- 

def readingFrameChecks ( proteins ):
    """checks conditions of start/stop amino acids and generates output for reading frames"""
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
        
# ---------------------------------------- 

# ---------------------------------------- 