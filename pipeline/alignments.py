from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist


def subsequence_without_gaps(query_seq, target_seq):
    return pairwise2.align.globalds(query_seq, target_seq, matlist.blosum62,
                                    -100, -100,
                                    penalize_end_gaps=False,
                                    one_alignment_only=True)


def align_possibly_gapped_sequence_on_its_complete_version(query_seq,
                                                           target_seq):
    return pairwise2.align.localxd(query_seq, target_seq, -1, 0, -100, -100,
                                   penalize_end_gaps=False,
                                   one_alignment_only=True)

def calc_identity(query_seq, target_seq):
  longer_seq = {}
  shorter_seq = {}

  if len(query_seq) >= len(target_seq):
    longer_seq = query_seq
    shorter_seq = target_seq
  else:
    longer_seq = target_seq
    shorter_seq = query_seq

  alignment_list = subsequence_without_gaps(longer_seq, shorter_seq)

  if not alignment_list:
      return 0

  alignment = alignment_list[0]

  longer_alignment = alignment[0]
  shorter_alignment = alignment[1]
  
  # calculate mismatches count
  mismatches_count = 0
  for i in range(len(longer_alignment)):
    if shorter_alignment[i] != longer_alignment[i]:
      mismatches_count += 1

  # count start and end gaps in the shorter alignment
  start_end_gaps_count = 0
  for i in range(len(longer_alignment)):
    if shorter_alignment[i] == '-':
      start_end_gaps_count += 1
    elif shorter_alignment[i] != '-':
      break
  for i in range(len(longer_alignment)):
    if shorter_alignment[len(longer_alignment) - i - 1] == '-':
      start_end_gaps_count += 1
    elif shorter_alignment[len(longer_alignment) - i - 1] != '-':
      break

  # account start and end gaps to calculate identity correctly
  mismatches_count -= start_end_gaps_count

  identity = (len(shorter_seq) - mismatches_count) / len(shorter_seq)

  return identity
