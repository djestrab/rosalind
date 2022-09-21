import numpy as np

from resources import RNA_CODON_TABLE


def load_string_from_file(filepath: str) -> str:
    with open(filepath) as f:
        lines = [line.strip() for line in f.readlines()]
        return "".join(lines)


def load_strings_from_file(filepath: str) -> list:
    with open(filepath) as f:
        return [line.strip() for line in f.readlines()]


def parse_fasta_from_file(filepath: str) -> dict:
    """
    Returns dictionary of DNAs parsed from FASTA format
    :param filepath:
    :return: dict of DNA strings
    """
    dna_dict = {}
    key = ""
    value = ""
    with open(filepath) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line[0] == ">":
                key = line[1:]
                value = ""
            else:
                value += line
            dna_dict[key] = value
    return dna_dict


def print_list(list_to_print: list):
    print(" ".join(map(str, list_to_print)))


def print_list_of_tuples(list_to_print: list):
    print("\n".join([str(i[0]) + " " + str(i[1]) for i in list_to_print]))


def print_profile_matrix(profile_matrix: dict):
    """
    Resulting format:
    A: 5 1 0 0 5 5 0 0
    C: 0 0 1 4 2 0 6 1
    G: 1 1 6 3 0 1 0 0
    T: 1 5 0 0 0 1 1 6
    """
    final_string = ""
    for key, value in profile_matrix.items():
        final_string += key + ": " + " ".join(map(str, value)) + "\n"
    print(final_string)


def count_nucleotides(dna: str):
    """
    ID: DMA
    Given: A DNA string s of length at most 1000 nt.
    Return: Four integers (separated by spaces) counting the respective number of times that the symbols
            'A', 'C', 'G', and 'T' occur in s.
    """
    return tuple(dna.count(nucleotide) for nucleotide in "ACGT")


def transcribe_dna(dna: str):
    """
    :param dna: A DNA string having length at most 1000 nt.
    :return: The transcribed RNA string of given DNA.
    """
    return dna.replace("T", "U")


def make_reverse_complement(dna: str):
    """
    :param dna: A DNA string s of length at most 1000 bp.
    :return: The reverse complement sc of s.
    """

    complements = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    return "".join([complements[n] for n in dna][::-1])


def count_fibonacci_rabbits(time: int, reproduction_number: int):
    """

    :param time: Positive integer n≤40
    :param reproduction_number: Positive integer ≤5
    :return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each
    generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair)
    """

    adult_pairs, young_pairs = 0, 1

    # every month:
    while time > 1:
        born_pairs = adult_pairs * reproduction_number
        adult_pairs += young_pairs
        young_pairs = born_pairs
        time -= 1

    return adult_pairs + young_pairs


def get_cg_content(dna) -> float:
    """

    :param dna:
    :return:
    """
    final_count = 0
    for char in dna:
        if char in "CG":
            final_count += 1
    return final_count / len(dna)


def get_highest_cg_content(dna_dict):
    """
    Return percentage of CG in DNA (0-100)
    :param dna_dict:
    :return:
    """
    highest_cg_content = 0
    highest_id = ""
    for dna_id, dna_string in dna_dict.items():
        cg_content = get_cg_content(dna_string)
        if cg_content > highest_cg_content:
            highest_cg_content = cg_content
            highest_id = dna_id
    return highest_id, round(highest_cg_content * 100, 5)


def calculate_hamming_distance(dna_1: str, dna_2: str) -> int:
    dna_length = len(dna_1)
    if len(dna_2) != dna_length:
        print("DNAs are not of the same length!")
        return -1

    distance = 0
    for idx in range(dna_length):
        if dna_1[idx] != dna_2[idx]:
            distance += 1
    return distance


def calculate_prob_dominant_allele(k: int, m: int, n: int) -> float:
    """
    k, m, n = population of organisms with individuals possessing different alleles for a factor:
    :param k: Homozygous dominant
    :param m: Heterozygous
    :param n: Homozygous recessive

    :return: The probability that two randomly selected mating organisms will produce an individual
    possessing a dominant allele (and thus displaying the dominant phenotype).

    Algorithm:
    1) Calculate probability of selecting individuals of given populations.
    2) Calculate probability that they wil produce offspring with dominant allele
    3) Multiply them
    4) Resulting probability is sum of the multiplied options
     
    """
    whole = k + m + n
    populations_probs = {   # tuple(prob of selecting individuals, prob of individuals having dominant offspring)
        "kk": (((k/whole)*((k-1)/(whole-1))), 1),
        "mm": (((m/whole)*((m-1)/(whole-1))), 0.75),
        "nn": (((n/whole)*((n-1)/(whole-1))), 0),
        "km": ((k/whole)*(m/(whole-1)), 1),
        "mk": ((m/whole)*(k/(whole-1)), 1),
        "kn": ((k/whole)*(n/(whole-1)), 1),
        "nk": ((n/whole)*(k/(whole-1)), 1),
        "mn": ((m/whole)*(n/(whole-1)), 0.5),
        "nm": ((n/whole)*(m/(whole-1)), 0.5)
    }
    prob = 0
    for probs_tuple in populations_probs.values():
        prob += probs_tuple[0] * probs_tuple[1]
    return round(prob, 5)


def translate_rna_into_protein(rna: str) -> str:
    """

    :rtype: object
    """
    codons = [rna[i:i+3] for i in range(0, len(rna), 3)]
    amino_acids = [RNA_CODON_TABLE[codon] for codon in codons if RNA_CODON_TABLE[codon] != "Stop"]
    protein = "".join(amino_acids)
    return protein


def find_motifs_locations(dna: str, motif: str) -> list:
    """

    :param dna:
    :param motif:
    :return: All locations (indexed from 1) of motifs as a substring of DNA.
    """
    return [idx+1 for idx in range(len(dna)) if motif == dna[idx:idx+len(motif)]]


def find_common_ancestor(dnas_pasta: dict) -> (str, dict):
    """
    Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) // in FASTA format.
    Return: A (one of) consensus string(s) and profile matrix for the collection.

    Alg:
    1) Make Matrix of DNAS
    2) Transpose
    3) Count char occurrences and write to profile matrix
    4) Create consensus string

    :param dnas_pasta:
    :return:
    """

    nucleotides = "ACGT"

    dnas_matrix = [list(dna) for dna in dnas_pasta.values()]
    trans_matrix = np.transpose(dnas_matrix)

    profile_matrix = {k: [0 for idx in range(len(dnas_matrix[0]))] for k in nucleotides}
    for char in nucleotides:
        for row_idx, row in enumerate(trans_matrix.tolist()):
            profile_matrix[char][row_idx] = row.count(char)

    consensus_string = ""
    for idx in range(len(dnas_matrix[0])):
        highest = ["", 0]
        for char, counts in profile_matrix.items():
            highest = [char, counts[idx]] if counts[idx] > highest[1] else highest
        consensus_string += highest[0]

    return consensus_string, profile_matrix


def count_mortal_fibonacci_rabbits(time: int, max_age: int) -> int:
    """
    ID: FIBD
    :param time: Positive integer n≤100
    :param max_age: Positive integer ≤5
    :return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months
    """

    ages = [1] + [0]*(max_age-1)    # list of rabbit pairs counts for every month
    while time > 1:
        births = sum(ages[1:])
        for age in range(max_age-1, 0, -1):
            if age == max_age-1:
                ages[age] = 0
            ages[age] = ages[age-1]
            ages[age-1] = 0
        ages[0] = births
        time -= 1

    return sum(ages)


def make_adjacency_list_naive(dnas: dict, overlap: int) -> list:
    """
    :param dnas: A collection of DNA strings in FASTA format having total length at most 10 kbp
    :return: The adjacency list corresponding to O-3. You may return edges in any order.
    """
    adj_list = []
    for name, dna in dnas.items():
        for name2, dna2 in dnas.items():
            if dna2[-overlap:] == dna[:overlap] and name2 != name:
                adj_list.append((name2, name))

    return adj_list


def calculate_expected_offspring(couples_distr: tuple) -> float:
    """
    Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the number
    of couples in a population possessing each genotype pairing for a given factor. In order, the six given integers
    represent the number of couples having the following genotypes:
    1: AA-AA
    2: AA-Aa
    3: AA-aa
    4: Aa-Aa
    5: Aa-aa
    6: aa-aa
    :param couples_distr: Six nonnegative integers
    :return: The expected number of offspring displaying the dominant phenotype in the next generation,
    under the assumption that every couple has exactly two offspring
    """
    assert len(couples_distr) == 6
    gen_probs = (1.0, 1.0, 1.0, 0.75, 0.5, 0)
    num_offspring = 2
    return sum([num_couples * num_offspring * prob for num_couples, prob in zip(couples_distr, gen_probs)])


def find_shared_motif(dnas: dict) -> str:
    """
    :param dnas: A collection of k (k≤100) DNA strings of length at most 1 kbp each in FASTA format.
    :return: A longest common substring (one of if more than one) of the collection.
    """

    def generate_substrings(input_string: str) -> str:
        """
        Yield substrings starting from the longest
        """
        for substring_length in range(len(input_string), 0, -1):
            start_idx = 0
            while start_idx + substring_length < len(input_string):
                yield input_string[start_idx:start_idx+substring_length+1]
                start_idx += 1

    # Prepare DNAs and sort them in ascending order
    dnas = list(dnas.values())
    dnas.sort(key=len)

    # Make substrings for the shortest DNA and iterate
    for substring in generate_substrings(dnas[0]):
        common_substring_found = False
        # Start comparing DNAs with substrings starting from the longest substring
        for compared_dna in dnas[1:]:
            if substring not in compared_dna:
                break
            common_substring_found = True
        # return first found substring found in all DNAs
        if common_substring_found:
            return substring
    return ""

