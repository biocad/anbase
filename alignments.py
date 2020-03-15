from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist


def subsequence_without_gaps(query_seq, target_seq):
    return pairwise2.align.globalds(query_seq, target_seq, matlist.blosum62,
                                    -100, -100,
                                    penalize_end_gaps=False,
                                    one_alignment_only=True)


def align_possibly_gapped_sequence_on_its_complete_version(query_seq,
                                                           target_seq):
    return pairwise2.align.localxs(query_seq, target_seq, -1, 0,
                                   penalize_end_gaps=False,
                                   one_alignment_only=True)
