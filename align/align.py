# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
       
        self._backtrack = None
        # # Init alignment_score
        self.alignment_score = 0

        # # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        self.MATCH = 0      # Diagonal move
        self.GAP_A = 1      # Vertical move (gap in sequence A)
        self.GAP_B = 2      # Horizontal move (gap in sequence B)



        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once

        self._seqA = seqA
        self._seqB = seqB

        m, n = len(seqA), len(seqB)
        
        # Initialize matrices
        # score_matrix[k][i][j] where k represents:
        # k=0: match/mismatch scores
        # k=1: gap in sequence A
        # k=2: gap in sequence B

        score_matrix = [[[float('-inf') for j in range(n+1)] for i in range(m+1)] for k in range(3)]
        score_matrix = np.array(score_matrix)  
        # inf for all

        # Backtrack matrix stores which matrix (0,1,2) gave the best score
        backtrack = [[0 for j in range(n+1)] for i in range(m+1)] # just 2D\
        backtrack = np.array(backtrack)
        
        # Initialize first cell
        score_matrix[0][0][0] = 0

        ### Initializations from video

        for i in range(1, m+1):
            score_matrix[1][i][0] = self.gap_open + (i-1) * self.gap_extend  # Gap in seqB
            #backtrack[i][0] = 1  # Indicates a gap in seqB (Insertion in seqA)

        for j in range(1, n+1):
            score_matrix[2][0][j] = self.gap_open + (j-1) * self.gap_extend  # Gap in seqA
            #backtrack[0][j] = 2  # Indicates a gap in seqA (Insertion in seqB)
            

        #print(f"Initial matrices\n M : \n {score_matrix[0]}\n Ix : \n {score_matrix[1]}\n Iy : \n {score_matrix[2]}")
        # Fill the matrices
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Calculate match/mismatch score: diagonal + substitution score
                match_score = score_matrix[0][i-1][j-1] + self.sub_dict[(seqA[i-1], seqB[j-1])]
                sub_score = self.sub_dict[(seqA[i-1], seqB[j-1])]
                
                # Calculate gap scores for sequence A
                # Open gap in A or extend gap in A
                gapA_extend = score_matrix[1][i-1][j] + self.gap_extend # pull from gap matrix (extending)
                gapA_open = score_matrix[0][i-1][j] + self.gap_open # pull from mismatch matrix
                score_matrix[1][i][j] = max(gapA_extend, gapA_open) # store gapped score in gap matrix A
                
                # Calculate gap scores for sequence B
                gapB_extend = score_matrix[2][i][j-1] + self.gap_extend
                gapB_open = score_matrix[0][i][j-1] + self.gap_open
                score_matrix[2][i][j] = max(gapB_extend, gapB_open)
                
                # Find the maximum score and store which matrix it came from
                scores = [match_score , score_matrix[1][i][j], score_matrix[2][i][j]]
 
                score_matrix[0][i][j] = max(scores) # store max score in match/mismatch matrix 
                #(we do this to pull the final total score)
                backtrack[i][j] = np.argmax(scores)

        self._backtrack = backtrack
        self._align_matrix = score_matrix
        #print("Finished populating matrices")
        #print(f"Populated matrices\n M : \n {score_matrix[0]}\n Ix : \n {score_matrix[1]}\n Iy : \n {score_matrix[2]}")
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        i, j = len(self._seqA), len(self._seqB)
        aligned_seqA, aligned_seqB = [], []
        score = self._align_matrix[0][i][j]
        
        while i > 0 or j > 0:
            if i > 0 and j > 0 and self._backtrack[i][j] == 0:
                # Match/Mismatch
                aligned_seqA.append(self._seqA[i - 1])
                aligned_seqB.append(self._seqB[j - 1])
                i -= 1
                j -= 1
            elif i > 0 and self._backtrack[i][j] == 1:
                # Gap in seqB (Insertion in seqA)
                aligned_seqA.append(self._seqA[i - 1])
                aligned_seqB.append('-')
                i -= 1
            elif j > 0 and self._backtrack[i][j] == 2:
                # Gap in seqA (Insertion in seqB)
                aligned_seqA.append('-')
                aligned_seqB.append(self._seqB[j - 1])
                j -= 1
            elif self._backtrack[i][j] == 1:
                aligned_seqA.append(self._seqA[i-1])
                aligned_seqB.append('-')
                i -= 1
            elif self._backtrack[i][j] == 2:
                aligned_seqA.append('-')
                aligned_seqB.append(self._seqB[j-1])
                j -= 1
        
        self.seqA_align = ''.join(reversed(aligned_seqA))
        self.seqB_align = ''.join(reversed(aligned_seqB))  
        self.alignment_score = score
        # Reverse sequences as we built them backwards
        return score, ''.join(reversed(aligned_seqA)), ''.join(reversed(aligned_seqB))



def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
