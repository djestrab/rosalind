from main import *
from pprint import pprint

# DNA
# print(count_nucleotides("GTGGCTTTTGAAGGCATCGTAGTACTCTTCCTAAAAAACTCCATCTTAAGACCTCCTCCCCTAGATACTAGAGTCACTTGCCTGCAGAGGATACTCGCTAATCTGTCATGATGTCGCGTGTGTCAAAGGATGTTCGCTGCGGACCAAGCCATCTCGGCGGTGCACTCATCGTCCTGCCTGTAGCGTGCGTTGGTGCCCCCTAACACTCGAGGGTTATGAGTAAGGTGGAAGAGACCTAACGAGAAGTTCCCTTCATCATCGTGAGATACTTGAATCCGTTGCACACGCGCGAACCATTTATACGATACTACTATACATGGCGAATTCAAGCCTAGCGTCGGACTGTTAATCGTTCTAGGGAAGGGCAGGGCCCTATCATAGCACAACTCCGCTTGGCGCGTCCTCGACGCTGGAAATGCGACCATGCGGGCGTTAGCGGGTAGATTTCCTTAATGACTTCATTTGCGCAGTTACACCCGGATCAATTCGCACACTGGTGCGCTCCACAGCCGTCGACTACAGGGGACGCTCACCCTAATCCCTCAAGCTCGGTTATGTGTCATGGAGTCCCGATTACGCATACCCCTAGCGGATTGCGCAAGCGTACCGCCGGACAAGAACCTTACCAGCCTACCAACGGTGGGTTGGAATGATCCAAGCAACTTGCCCGCCTGGGCCATTTAGGACGAATACTATGTGAATCAGTCACGAGTGGACAGTGCTATGAGACCTGATTAATACATACGCATTTTTTCTGACACTGGCTCCTGACATCGGCCCTCATGCCTTTAAAAGCTATTAGTCCGTGCACAGCACGTCGAGAGTTTTAGTGCAGCATCCCCGGTTCCCGGAGCAACTTTTGCTGCTTTAGGGCCTACACGGG"))

# RNA
# print(transcribe_dna("CTTTCTTAGAATCTGGATGACTAGAGTTACTACAACAACGCTCAACGTATATCTGCTGGAGGCGAACAGGCGGCGGATGGCTGCACCACCGCATCAGCCCAGTAGTGTGCCAGACCAGGCTCTGTAGGCTTCTCACTTCGACTCGCATGCTGTCCCCTTAAGCACGGCGCATTGGGAGAGATAGCGCCCCGGACGTCACCCCGGCATTCATAAAGCAATAGTGACGTACGAAGTACCCGAGAGGAGGGAAAGGCCCTCCTATTACGTGCGACAGCACTGCGTATCTAGAGCTCCGACGTTGTCAGCCCAGGATATAGTTCACGCGTCGCGAACACCGAGGGTTAGTATATACAATGTCCCTCGCTCTTCGGAGCGATTACATCATCACTCGAGGTACCATAGAATACTAGTCCACAGATCAAGAAGGTAAGTGACGACGGGTAGCAGGCCCCCGCCTGCACCGAAAGCGGTGCCGTAAACCCAGCGCACAGATGAGATCATATGCAATGCGAGCGGAGCCACAATCTTACTCGGACCACCTCTAGTCCAACGGGATGGGGTAGCAAAAGTCATGCCCCAGGTTCTTTTCGGCTCCACCAGCCGACAAGTGCGCGCCCTATCTCCCAGACCGGTGGATCTGGCCGACTGATTACTCCTTCCTCTGGCACGACGTTGAACCAATAGCATGAGCTGAATGTTCCAGTCCGTTTAACCAGCCCAATGATACTATGCGGTACAGAGCGGACATTTATAGGATGGAACCCGTTATGCCACGATGTTCGCACACTACGGGAGTCATCGGATTTCCAGTCGTTGCCACAATCAGCTTTATGGATCATGTTCCGGGGCATCGGGCTGACTCTAGGCCAGCGGCGTACTCAAAGACAGACAATGTAGAGCTAACTTAAGTTTCAACTACGCACGAATCCGCAACCATT"))

# REVC
# print(make_reverse_complement("AGAAAATCCCGCTCATGCAAGGCGATTTAGCTTTTTGATAACACTGATCTTCACGTATAGAACGCTAACCCTAGTGGCCCGCACACGGGCGTCGCAATCAGGAGTACCCTGATCATAGTAAGGTCACACGAGCCCCACCAAACATCATCGCTTATTGCAATAGAATAGATTACCATAGAGCGTCCACCGGGTGCTTCTGCACATCTCGTATGCAAGTCCCTCGTTGTCTTAATATTACTGAGCCGGTGAACCAGGTGTCCGTCAGTCGCAGGTTGATGGTTTTTTAGGGGCGCTTGGAACCACACCGTCAATTACTGCGCCAGGAATGCCTGCACCAATCCTTAACCCGAATCATCCTATTTAGTTTTTGAGGAAACAGGGTTTGGAAGATCTGGTGAGATCTTTGCATCAAGAATTCAGAAGCGTGCGAATTACGAGCCAACCTACTTCATGAATGTTGAAACTATTCATGGTCCGATGGCCTCCGGCTCCAAGGGTTACACTAGGACACGACGAGGCACAAATAGTACGGTAGGCTGGATCGTAGATGCACAGCGAGGGGTATCTACCTATGCGGCCACTCAAAGGGACAGCGTATTATGATCGCCCTAACATACATATGTTGGGTAACTCAGATCAGAGTCTGCGGCCACGCATCGAGCGCCTGATGACCGTCGGCGCTGGAGGGGGGAGGAGCTTAGGTACACCATCATCGTACATGCAAGTTATACGGGCTAGGGATAACTTCCACCTTTCGCAACTAACATCTGTGAGTTCATCCTTTCCAACTCTCTCATCGCGGATGCTTGCGTATGAACGAAAATTTGGGATAACATATAACTCCAGATTGACTCCTTATCCCCCGACGCAGTTAATACTGAAGGTGTGAAAATACCCTATACCTTAGCGTCACCTTAGGCAGGTGCGTCTTA"))

# FIB
# print(count_fibonacci_rabbits(36, 2))

# GC
# print(get_highest_cg_content(parse_pasta_from_file("C:/Users/danie/Downloads/rosalind_gc.txt")))

# HAMM
# print(calculate_hamming_distance("CATGAGGCTTAAAGATCCGAGCGAGGCTAGAACACCTTCCGTCATTTCGCATGAGCCTGATATCAATTGTCCCTTACAAGCCACTGTAGAGCGGAGGCATATTGCTCTTGAAAAGCATGAGCAACGTACATTAACACAACTTTCCACACGAACAGTATGACACCCGGTACCCATAATAACCATTATCAACCTGAATGTGACCGTCATACCACTCTACCCACAGCAGTTTAGGGTGAGGCTTTACATCGGCGGAGTGAGGACCCCCATGACCAAAGGAGCAGGAACCCGCCATATCACGTTGAAATATTGGCGTGATGCTCAAATTGTTTTGTACGAGTCCCATTTAACAAGCGGTACCACCTGGTAGGTCTTTGACCGTTTGCGGCCTCTGAATCCGTTCGGCAGAAAATAACGACATCTAGCCTCAGGGTGTGGACTCGTCATAAGTCTGCCGATAACCTACCTAGTGCAAAGGAGAACCTTACGAACCTGAGGTTCCCTTCTTGGTTCAGATTTATGTTACCTCACCGTCAGAGCGGGGCGGACGCCTGATTAGACAGTTTGGCTGCCAATCCGCTTGACTTTTCACCTTCCTAAACTCGCACGCCACCGCGGTAAACCGCCTTAGAGGCATAACTCAAGACATAGGGCAATCGTAGGTCGTTGGGGGGAAAGTCCAAGAGAGTTAGAGAGACGCTCCCATCCGTTGTTCAGAAGACCGATCAAAAATCTGCAAGCGTACTTTTTTCAGTCTCTGGGTATGGTACGCTTATGCGGAGTAGCAGGCCAAGCCTGGAGCAGCTACGCAGAAAGTATATGAGTAGCTGAAAAATCGTAAGCATGACCCACAAACCTGACATTGGACAATACACTAAATTTCGGCCTCGGTAGCAAGCAAAGAGGGTGTTA", "CTTCGAGCACACGATAACAAGCGCGTTTTACAAATTATCTGCTGTTTAGAGTGCCCCAGCAAGCTCATTTACTTTGCGATACCGTCTATAGCGGAGGTAATTTGATCAATAAACTTATCTGCAACCGCTCATAAGAACCCCTCTCACACGCTCAGAAGGACATCGGGTAACCCGATAGATTATTGAACACACGAGTGAATTATATATTTGAGTTCAATTTCACCTCTTTTGCGTGATGCACTAACGCGGCCCAGTGGGTGCGAGCTTTACCGATCCTGGCGCAACTGGCTTTCCCAAGTTTAAATACCGGCCGTTAACGCAAATGCGCTCGTAGTAGTCCGAGTTACTAACGGGTACCCGCTATAAGGCCTTACCCGGTTTTCGGATTGTCTACCCTGACGGAGGATCTTAACTACATATAAGCTCAGCAATTGGTCCGGCACTAAGTCCACTTCGGACTTTCATACTCGCACGTGGAACTATTCAAGCTCAATCTTTCAATGTGAGGTCGTGCACCTAATACCTCGCCTTGAACAGTGCGGCGAGGACTGATTAGCCAGGTTAGCCGAGAACTAACTCGACGCTATGGCTTTCTACCGTAACACGCCAGCGGGCTTGTGCACATCCTAGTGGGTTCTAACGATGGAGCGCCTGCCTTTGTCTTTTTGGGTTGAGTGCTGCAGTGTCATAGGGCCGGATCCGCTCTTTATTGCGAAGACGGTAAGAAGCTTTGCAACGTCTAGTCGTCCACTCTCGTGAGGTCAGAGGCGATGACGGAGTGGCAGTCCTTGCCCTTGTCCACACATCGAAGGGCATAATTCTTAATACTTAATCAATTGCGTCACCGTCCAAGCTATGACCCCACAGTACACTATACTTATTAGTAGCTAGGATACACCCCGCACTATA"))

# PRB
# print(calculate_prob_dominant_allele(15, 30, 16))

# PROT
# print(translate_rna_into_protein(load_string_from_file("C:/Users/danie/Downloads/rosalind_prot.txt")))

# SUBS
# rosalind_file = "C:/Users/danie/Downloads/rosalind_{}.txt".format(task_id.lower())
# print_list(find_motifs_locations(*load_strings_from_file(rosalind_file)))

# CONS
# rosalind_file = "C:/Users/danie/Downloads/rosalind_{}.txt".format(task_id.lower())
# ble = find_common_ancestor(parse_fasta_from_file(rosalind_file))
# print(ble[0])
# print_profile_matrix(ble[1])

# FIBD
# print(count_mortal_fibonacci_rabbits(99, 19))

# GRPH
# rosalind_file = "C:/Users/danie/Downloads/rosalind_{}.txt".format(task_id.lower())
# print_list_of_tuples(make_adjacency_list_naive(parse_fasta_from_file(rosalind_file), 3))

# IEV
# print(calculate_expected_offspring((17152, 19940, 16855, 17348, 17629, 19269)))

# LCSM
# rosalind_file = "C:/Users/danie/Downloads/rosalind_{}.txt".format("LCSM")
# print(find_shared_motif(parse_fasta_from_file(rosalind_file)))







