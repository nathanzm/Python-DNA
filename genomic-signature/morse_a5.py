"""----------------------------------------------------------------------------
SUMMARY: This program generates a report that serves as a preliminary study to 
compare and contrast certain features of sequences between any given amount of
different organisms, or just one organism. By using the analysis of different
sizes of motifs within a DNA sample, this program can give the user an idea of 
an organism's genomic signature, which helps the user determine characteristic
frequencies of oligonucleotides in the given genome or sequence.

Programmer: Nathan Morse
    
INPUT: There are two sources of input for this program:
1. User input for the L-mer size (in base pairs) that the user wants to use for
the size of motifs that will be analyzed using the 2nd input mentioned below.
This input should be an integer between 4 and 8 (inclusively) and if it is not, 
the user will be prompted to enter an input again. The minimum and maximum size
of the L-mer that can be entered can be changed by changing the variables that
are in the main at the bottom of the program.
2. The second input is a file directory containing input files for analysis.
This file directory should contain only FASTA formatted files of entire genomes
or just individual genes. Any files that are not FASTA formatted in the directory
will simply be ignored. The directory can contain as many files as the user wants.
The directory name can be changed by changing the variable "inputFilesDirectory"
in the main at the bottom of the program.

OUTPUT: There are two forms of output for this program:
1. The 1st output the user will notice is directly to the console. A title of
this report and an app name for the program with an ASCII design along with a
prompt for the user to enter an integer for the L-mer size between the given 
minimum and maximum L-mer size (this can be changed by changing the variables in
the main). If the user does not enter a correct integer, the user will be receive
an according error message. Once the user enters a correct input, the files within
the input directory are printed to the screen with their according total DNA 
lengths (in bp). If the input directory does not exist or does not contain any 
FASTA formatted files, an appropriate error message will be printed. The last
successful line printed to the screen is a message telling the user that an 
output file has been successfully generated.
2. The 2nd output is an excel file that is written to the same directory as this
program. This comma separated value file (.csv) contains each input file's name,
total number of motifs based on the size of L-mer entered by the user, the top
number of motifs (determined by the variable "numberOfTopRankings" in the main;
the default is the top 10) with their according name, frequency, and proportion.
-------------------------------------------------------------------------------
"""
import BioDNA

#--------------------------------------------------------------------------
def PrintProgramTitle():
    # this function prints the title of the program to the screen
    print("==================================================================")
    print("=-=-=-=-=-=-=-=-= + Genomic Signature Analyzer + =-=-=-=-=-=-=-=-=")
    print("    O       o O       o O       o O       o O       o O       o")
    print("    | o   O | | o   O | | o   O | | O   o | | O   o | | O   o |")
    print("    | | O | | | | O | | | | O | | | | O | | | | O | | | | O | |")
    print("    | O   o | | O   o | | O   o | | o   O | | o   O | | o   O |")
    print("    o       O o       O o       O o       O o       O o       O")
    print("            GSA (build v1.3.12) - by Nathan Morse")
    print("==================================================================")
#-----------------(end of PrintProgramTitle())-----------------------------

def main():
    MIN_LmerSize = 4 # minimum Lmer size
    MAX_LmerSize = 8 # maximum Lmer size
    inputFilesDirectory = "inputDirectory" # directory of input files
    numberOfTopRankings = 10 # number of top rankings to display
    # variable prompt is for printing correct statement for user input based on LmerSize variables:
    prompt = "\nEnter an L-mer size between "+str(MIN_LmerSize)+" and "+str(MAX_LmerSize)+" (inclusively): "
    
    PrintProgramTitle() # call to function that prints title of program
    
    """User input size of Lmer and check if it is an integer between the
    MIN_LmerSize and MAX_LmerSize inclusively:
    """
    while (True): # keeping going through while loop until user input is valid
        try: # if the user doesn't enter an integer type, immediately go to ValueError
            LmerSize = int(input(prompt))
            # check if LmerSize between min and max and if not, ask for user input again:
            while (MIN_LmerSize > LmerSize or MAX_LmerSize < LmerSize):
                print("You did not enter a valid L-mer length, try again!")
                LmerSize = int(input(prompt))            
        except ValueError: # print error if input is not an integer
            print("The input you have entered is not an integer!")
            continue # go back to try above
        else: # break out of while loop once input is valid
            break
    
    print("You successfully entered", LmerSize, "for the size of the L-mer.\n")
    
    # call to inputDirectoryLengths on argument that is the file path:
    BioDNA.inputDirectoryLengths(inputFilesDirectory)
    
    """ call to breakIntoMotifs on arguments that include the file directory, L-mer size,
    and the number of top rankings that you want to be written to the output file:
    """
    BioDNA.breakIntoMotifs(inputFilesDirectory, LmerSize, numberOfTopRankings)
     
# --- end main() ---------

#---------------------------------------------------------
# Python starts here ("call" the main() function at start
if __name__ == '__main__':
    main()
#---------------------------------------------------------