import pytest
from main import *


def test_print_list():
    assert print_list([2, 4, 5, 52222, "fefw"]) == "2 4 5 52222 fefw"


def test_count_nucleotides():
    assert count_nucleotides("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC") == (20, 12, 17, 21)


def test_transcribe_dna():
    assert transcribe_dna("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU"


def test_make_reverse_complement():
    assert make_reverse_complement("AAAACCCGGT") == "ACCGGGTTTT"


def test_count_fibonacci_rabbits():
    assert count_fibonacci_rabbits(5, 3) == 19


def test_parse_pasta_from_file():
    pasta_dict = {
        "Rosalind_6404": "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
        "Rosalind_5959": "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
        "Rosalind_0808": "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
    }
    assert parse_fasta_from_file("files/fasta.txt") == pasta_dict


def test_get_cg_content():
    assert get_cg_content("AAAACCCCGGGGTTTT") == 0.5


def test_get_highest_cg_content():
    # TODO: use fixtures
    fasta_dict = {
        "Rosalind_6404": "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG",
        "Rosalind_5959": "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC",
        "Rosalind_0808": "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
    }
    assert get_highest_cg_content(fasta_dict) == ("Rosalind_0808", 60.919540)


def test_calculate_hamming_distance():
    dna_1 = "GAGCCTACTAACGGGAT"
    dna_2 = "CATCGTAATGACGGCCT"
    assert calculate_hamming_distance(dna_1, dna_2) == 7


def test_calculate_prob_dominant_allele():
    assert calculate_prob_dominant_allele(2, 2, 2) == 0.78333


def test_translate_rna_into_protein():
    assert translate_rna_into_protein("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA") == "MAMAPRTEINSTRING"


def test_find_motifs_locations():
    assert find_motifs_locations("GATATATGCATATACTT", "ATAT") == [2, 4, 10]


def test_find_common_ancestor():
    profile_matrix = {
        "A": [5, 1, 0, 0, 5, 5, 0, 0],
        "C": [0, 0, 1, 4, 2, 0, 6, 1],
        "G": [1, 1, 6, 3, 0, 1, 0, 0],
        "T": [1, 5, 0, 0, 0, 1, 1, 6]
    }
    assert find_common_ancestor(parse_fasta_from_file("files/fasta_cons.txt")) == ("ATGCAACT", profile_matrix)


def test_count_mortal_fibonacci_rabbits():
    assert count_mortal_fibonacci_rabbits(6, 3) == 4


def test_make_adjacency_list_naive():
    expected = [
        ("Rosalind_0498", "Rosalind_2391"),
        ("Rosalind_0498", "Rosalind_0442"),
        ("Rosalind_2391", "Rosalind_2323")
    ]
    assert frozenset(make_adjacency_list_naive(parse_fasta_from_file("files/fasta_grph"), 3)) == frozenset(expected)


def test_calculate_expected_offspring():
    assert calculate_expected_offspring((1, 0, 0, 1, 0, 1)) == 3.5


def test_find_shared_motif():
    assert find_shared_motif(parse_fasta_from_file("files/fasta_lcsm.txt")) == "?"

