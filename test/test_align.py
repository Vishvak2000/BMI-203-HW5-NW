# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    nw = NeedlemanWunsch(gap_open=-10, gap_extend=-1,sub_matrix_file="./substitution_matrices/BLOSUM62.mat")
    nw.align(seq1, seq2)
    assert np.allclose(nw._align_matrix[0],np.array(([[  0., -np.inf, -np.inf, -np.inf],
                                                    [-np.inf,   5.,  -5.,  -6.],
                                                    [-np.inf,  -5.,   4.,  -6.],
                                                    [-np.inf,  -6.,   0.,   5.],
                                                    [-np.inf,  -7.,  -5.,   5.]])))
    assert np.allclose(nw._align_matrix[1],np.array([[-np.inf, -np.inf, -np.inf, -np.inf],
                                                    [-10., -np.inf, -np.inf, -np.inf],
                                                    [-11.,  -5., -15., -16.],
                                                    [-12.,  -6.,  -6., -16.],
                                                    [-13.,  -7.,  -7.,  -5.]]))
    assert np.allclose(nw._align_matrix[2],np.array([[-np.inf, -10., -11., -12.],
                                                    [-np.inf, -np.inf,  -5.,  -6.],
                                                    [-np.inf, -np.inf, -15.,  -6.],
                                                    [-np.inf, -np.inf, -16., -10.],
                                                    [-np.inf, -np.inf, -17., -15.]]))
    
    

def test_nw_backtrace():
 
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch(gap_open=-10, gap_extend=-1,sub_matrix_file="./substitution_matrices/BLOSUM62.mat")
    nw.align(seq3, seq4)
    assert nw.seqA_align == "MAVHQLIRRP"
    assert nw.seqB_align == "M---QLIRHP"
    





