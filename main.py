# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    nw = NeedlemanWunsch(gap_open=-10,gap_extend=-1,sub_matrix_file="./substitution_matrices/BLOSUM62.mat")

    species = {
        "Gallus gallus": gg_seq,
        "Mus musculus": mm_seq,
        "Balaeniceps rex": br_seq,
        "Tursiops truncatus": tt_seq
    }
    #print("loaded in fastas")

    scores = {}
    for species_name, seq in species.items():
        score, _, _ = nw.align(hs_seq, seq)
        scores[species_name] = score
        #print(score)
    
    sorted_species = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    
    print("Species in order of most similar to human BRD2:")
    for species_name, score in sorted_species:
        print(f"{species_name}: {score}")

if __name__ == "__main__":
    main()
