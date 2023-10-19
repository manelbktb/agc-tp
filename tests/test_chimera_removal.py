"""Tests chimera removal"""
import pytest
import os
import nwalign3 as nw
from pathlib import Path
from .context import agc
from agc import get_chunks
from agc import get_unique
from agc import common
from agc import cut_kmer
from agc import get_unique_kmer
from agc import search_mates
from agc import detect_chimera
from agc import get_identity
from agc import chimera_removal


def test_get_chunks():
    """
    """
    # let's cut this 222nt sequence in chunks of 80 nt
    seq = "TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGGAAAGCGCA"
    try:
        chunks = get_chunks(seq, 80)
    except ValueError:
        # Congrats error correctly raised
        assert(True)
    chunks = get_chunks(seq, 50)
    # we should obtain 4 chunks of the correct size
    assert(len(chunks) == 4)
    # a 80nt chunk
    assert(chunks[0] == seq[0:50]) 
    # a 80nt chunk
    assert(chunks[1] == seq[50:100])

def test_unique():
    res = get_unique([1,2,3,2,3,4])
    assert(len(res) == 4)

def test_common():
    res = common([1,2,3], [2,3,4])
    assert(1 not in res)
    assert(2 in res)
    assert(3 in res)
    assert(4 not in res)


def test_cut_kmer():
    kmer_reader = cut_kmer("TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAG", 33)
    assert(next(kmer_reader) == "TGGGGAATATTGCACAATGGGCGCAAGCCTGAT")
    assert(next(kmer_reader) == "GGGGAATATTGCACAATGGGCGCAAGCCTGATG")
    assert(next(kmer_reader) == "GGGAATATTGCACAATGGGCGCAAGCCTGATGC")
    assert(next(kmer_reader) == "GGAATATTGCACAATGGGCGCAAGCCTGATGCA")
    assert(next(kmer_reader) == "GAATATTGCACAATGGGCGCAAGCCTGATGCAG")
    # next step should crash
    try:
        next(kmer_reader)
    except StopIteration:
        assert(True)


def test_get_unique_kmer():
    """
    """
    kmer_dict = get_unique_kmer({}, "TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAG", 0, 8)
    kmer_dict = get_unique_kmer(kmer_dict, "GGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGC", 1, 8)
    assert(len(kmer_dict) == 31)
    assert(len(kmer_dict["TGGGGAAT"]) == 1)
    assert(len(kmer_dict["GGGGAATA"]) == 2)
    assert(len(kmer_dict["GATGCAGC"]) == 1)


def test_search_mates():
    """
    """
    kmer_dict = get_unique_kmer({}, "TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAG", 0, 8)
    kmer_dict = get_unique_kmer(kmer_dict, "GGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGC", 1, 8)
    kmer_dict = get_unique_kmer(kmer_dict, "GGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCC", 2, 8)
    best_mates = search_mates(kmer_dict, "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCA", 8)
    assert(len(best_mates) == 3)
    assert(best_mates[0] == 2)
    assert(best_mates[1] == 1)
    assert(best_mates[2] == 0)


def test_detect_chimera():
    # best:
    # 1 , 1, 1, 0
    assert(not detect_chimera([[73.58, 73.79], [72.97, 77.06], [77.36, 80.58], [78.43, 78.43]]))
    assert(not detect_chimera([[62.6, 94.17], [62.6, 94.17], [62.6, 94.17], [62.6, 94.17]]))
    assert(detect_chimera([[98.0, 60.0], [100.0, 65.0], [100.0, 63.0], [64.0, 100.0]]))
    S000387216 = "GGAGGCTCGTACCGCTGTCTTGTTAAGGACTGGTTTTTTACTGTCTATACAGACTCTTCATACTACTGGATATCCTGATATGCGTTCGGATCGATTGTTGCCGTACGCTGTGTCGATTAAAGGTAATCATAAGGGCTTTCGACTTACGACTC"
    chimera_AJ007403 = "AAGACGCTTGGGTTTCACTCCTGCGCTTCGGCCGGGCCCGGCACTCGCCACAGTCTCGAGCGTCGTCTTGATGTTCACATTGCGTTCGGATCGATTGTTGCCGTACGCCTGTGTCATTAAAGGTAATCATAAGGGCTTTCGACTTACGACTC"
    S000001688 = "AAGACGCTTGGGTTTCACTCCTGCGCTTCGGCCGGGCCCGGCACTCGCCACAGTCTCGAGCGTCGTCTTGATGTTCACATGTAACGATCGCTTCCAACCCATCCGGTGCTGTGTCGCCGGGCACGGCTTGGGAATTAACTATTCCCAAGTCT"
    chunk_chim = get_chunks(chimera_AJ007403, 37)
    chunk_seq_list = [get_chunks(S000387216, 37)]
    chunk_seq_list += [get_chunks(S000001688, 37)]
    perc_identity_matrix = [[] for c in range(len(chunk_chim))]
    for i in range(len(chunk_seq_list)):
        for l,chunk in enumerate(chunk_chim):
            perc_identity_matrix[l].append(get_identity(
                        nw.global_align(chunk, chunk_seq_list[i][l], 
                            gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '../agc')) + "/MATCH")))
    # [[48.89, 100.0], [51.16, 100.0], [84.62, 58.54], [94.74, 45.65]]
    assert(detect_chimera(perc_identity_matrix))


def test_chimera_removal():
    chimerafree = chimera_removal(Path(__file__).parent / "test_sequences.fasta.gz",
        200, 3, 50, 8)
    assert(next(chimerafree)[0] == "ACTACGGGGCGCAGCAGTAGGGAATCTTCCGCAATGGACGAAAGTCTGACGGAGCAACGCCGCGTGTATGAAGAAGGTTTTCGGATCGTAAAGTACTGTTGTTAGAGAAGAACAAGGATAAGAGTAACTGCTTGTCCCTTGACGGTATCTAACCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAGTTAGTGGGCGTAAAGCGCGCGCAGGCGGTCTTTTAAGTCTGATGTCAAAGCCCCCGGCTTAACCGGGGAGGGTCATTGGAAACTGGAAGACTGGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGATATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGCGCGAAAGCGTGGGGAGCAAA")
    assert(next(chimerafree)[0] == "TAGGGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACATATGTGTAAGTAACTGTGCACATCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTACAGCGCG")

