CONTENTS:
- This README text file
- BioDNA.py : python file with defined useful functions
- morse_a3main.py : python file that executes functions in BioDNA.py

SUMMARY: By running the python file morse_a3main.py (assuming that BioDNA.py remains
in the same directory), the user can achieve three different reading frames with 
information from of an mRNA based off of a given DNA sequence and its anti-sense strand.
The information given gives the user an understanding of how the translation process
from mRNA to amino acids work and where/whether or not there are start and stop sites
within the protein (chain of amino acids). The user can easily determine by looking at the
output which reading frame is best.
    
Programmer: Nathan Morse
    
INPUT: A string of DNA (indicated by the variable "DNA" in morse_a3main.py)

OUTPUT: The DNA strand along with its complementary anti-sense strand and the transcribed
messenger RNA. In addition, there are three reading frames (which are achieved by shifting
the RNA over one character each time), all of which include the RNA length, the mRNA strand,
its aligned and translated protein, whether or not there is a Methionine amino acid at the
start site or elsewhere in the protein, and whether or not there is an amino acid signaling
a STOP at the end of the protein or elsewhere. 