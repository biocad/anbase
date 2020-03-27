import os
import string
from collections import defaultdict
from xml.etree import ElementTree

import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.StructureBuilder import StructureBuilder

from fetch_unbound_data import fetch_all_sequences, AG, AB, DB_PATH, DOT_PDB, \
    get_while_true, comp_name_to_pdb_and_chains, CHAINS_SEPARATOR, \
    is_subsequence_of, get_real_seqs, fetch_struct

FILTERED_STRUCTURES_CSV = 'filtered_for_unboundness.csv'
REJECTED_STRUCTURES_CSV = 'rejected_for_unboundness.csv'

FILTERED_COMPLEXES_CSV = 'filtered_complexes.csv'
REJECTED_COMPLEXES_CSV = 'rejected_complexes.csv'


def process_csv(csv):
    data = defaultdict(list)

    for i in range(len(csv)):
        data[csv.iloc[i]['comp_name']].append((csv.iloc[i]['candidate_type'],
                                               csv.iloc[i]['candidate_pdb_id'],
                                               csv.iloc[i][
                                                   'candidate_chain_ids']))

    return data


def fetch_number_of_assemblies(pdb_id):
    curl = 'https://www.rcsb.org/pdb/rest/bioassembly/' \
           'nrbioassemblies?structureId={}' \
        .format(pdb_id)

    r = get_while_true(curl)
    xml = ElementTree.fromstring(r)

    if 'count' not in xml.attrib:
        return fetch_number_of_assemblies(pdb_id)

    return int(xml.attrib['count'])


def fetch_all_assemblies(pdb_id):
    n = fetch_number_of_assemblies(pdb_id)

    res = []

    for i in range(n):
        curl = 'https://files.rcsb.org/download/{}.pdb{}'. \
            format(pdb_id, str(i + 1))

        path_to_tmp = os.path.join(DB_PATH, pdb_id + '_BA_' + str(i) + DOT_PDB)

        if os.path.exists(path_to_tmp):
            res.append(path_to_tmp)
            continue

        r = get_while_true(curl)

        if r is None:
            continue

        with open(path_to_tmp, 'w') as f:
            f.write(r)

        res.append(path_to_tmp)

    return res


def union_models(struct):
    models = list(struct.get_models())

    if len(models) < 2:
        return struct

    sb = StructureBuilder()

    sb.init_structure('ba')
    sb.init_model(0)

    chain_names = {}

    for model in models:
        for chain in model:
            chain_id = chain.get_id()
            name = chain_id

            if chain_id not in chain_names:
                chain_names[chain_id] = 0
            else:
                chain_names[chain_id] += 1
                name = chain_id + '_' + str(chain_names[chain_id])

            chain_copied = chain.copy()
            chain_copied.id = name

            sb.model.add(chain_copied)

    return sb.structure


class AssemblyMatchInfo:
    def __init__(self, assembly_id, matching, reason_bad=None):
        self.is_good = reason_bad is None
        self.id = assembly_id
        self.matching = matching
        self.reason_bad = reason_bad

    def __repr__(self):
        return str((self.is_good, self.id, self.matching,
                    self.reason_bad))


def check_structure(source_pdb_id, source_chain_ids, target_pdb_id, type):
    is_ab = type == AB
    pdb_parser = PDBParser()

    try:
        struct = pdb_parser.get_structure(source_pdb_id,
                                          fetch_struct(source_pdb_id))
    except Exception as e:
        print('ERROR:', e, flush=True)
        return []

    source_seqs = list(filter(lambda x: x[0] in source_chain_ids,
                              fetch_all_sequences(source_pdb_id).items()))
    real_source_seqs = list(
        zip([x[0] for x in source_seqs], get_real_seqs(struct, source_seqs)))

    target_seqs = fetch_all_sequences(target_pdb_id)

    res = []

    n = -1
    for assembly_path in fetch_all_assemblies(target_pdb_id):
        n += 1

        try:
            assembly_structure = pdb_parser.get_structure('ba', assembly_path)
        except Exception as e:
            print('BAD ASSEMBLY FOR {}:'.format(target_pdb_id), e, flush=True)
            continue

        assembly = union_models(assembly_structure)

        chains_in_assembly = [x.get_id().split('_')[0]
                              for x in assembly.get_chains()]

        assembly_ids_seqs = list(
            map(lambda x: (x, target_seqs[x]), chains_in_assembly))

        chain_matching = defaultdict(list)

        for chain_id, chain_seq in real_source_seqs:
            for target_chain_id, target_seq in assembly_ids_seqs:
                if is_subsequence_of(chain_seq, target_seq, is_ab=is_ab):
                    chain_matching[chain_id].append(target_chain_id)

        lens_of_matches = list(map(lambda x: len(chain_matching[x]),
                                   source_chain_ids))

        if not lens_of_matches:
            continue

        n_plus_one = n + 1

        if all(map(lambda x: x == 1, lens_of_matches)) and len(
                assembly_ids_seqs) == len(source_seqs):
            # assembly contains only matching with needed seqs

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching))
        elif all(map(lambda x: x == 1, lens_of_matches)) and len(
                assembly_ids_seqs) != len(source_seqs):
            # assembly contains matching and some other chains

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching,
                                         reason_bad='additional_chains'))
        elif lens_of_matches[0] > 0 and all(
                map(lambda x: x == lens_of_matches[0],
                    lens_of_matches)) and len(assembly_ids_seqs) == \
                lens_of_matches[0] * len(source_seqs):
            # assembly contains potential homomer that contains many matchings

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching,
                                         reason_bad='potenial_homomer'))
        elif lens_of_matches[0] > 0 and all(
                map(lambda x: x == lens_of_matches[0],
                    lens_of_matches)) and len(assembly_ids_seqs) != \
                lens_of_matches[0] * len(source_seqs):
            # assembly contains potential complex homomer that
            # contains many matchings

            res.append(AssemblyMatchInfo(n_plus_one, chain_matching,
                                         reason_bad='potential_'
                                                    'complex_homomer'))

    return res


def get_pdb_ids(l, ty):
    return list(
        frozenset(map(lambda x: x[1], filter(lambda x: x[0] == ty, l))))


def matching_to_str(chains, matchings):
    n_matchings = len(matchings[list(matchings.keys())[0]])

    by_matching = []

    for i in range(n_matchings):
        by_matching.append([])
        for chain in chains:
            by_matching[-1].append(matchings[chain][i])

    return '|'.join(map(lambda x: CHAINS_SEPARATOR.join(x), by_matching))


def filter_candidates_pack(comp_name, pdb_id, candidate_pdb_ids, chain_ids, ty,
                           filtered_csv, rejected_csv):
    for candidate_pdb_id in candidate_pdb_ids:
        chains_str = CHAINS_SEPARATOR.join(chain_ids)

        assemblies = check_structure(pdb_id,
                                     chain_ids,
                                     candidate_pdb_id, ty)
        for assembly in assemblies:
            matching_str = matching_to_str(chain_ids, assembly.matching)

            if assembly.is_good:
                filtered_csv.write(','.join(
                    [comp_name, ty, chains_str, candidate_pdb_id,
                     matching_str, str(assembly.id)]) + '\n')
            else:
                rejected_csv.write(','.join(
                    [comp_name, ty, chains_str, candidate_pdb_id,
                     matching_str, str(assembly.id),
                     assembly.reason_bad]) + '\n')

        filtered_csv.flush()
        rejected_csv.flush()


def filter_for_unboundness(processed_csv):
    post_processed = set([])

    if os.path.exists('post_processed.csv'):
        with open('post_processed.csv', 'r') as f:
            for line in f.readlines():
                post_processed.add(line.strip())

    mode = 'a' if post_processed else 'w'

    with open(FILTERED_STRUCTURES_CSV, mode) as filtered_csv, open(
            REJECTED_STRUCTURES_CSV, mode) as rejected_csv, open(
        'post_processed.csv', 'a') as post_processed_csv:

        if mode == 'w':
            filtered_csv.write(
                'comp_name,type,chain_ids,candidate_pdb_id,'
                'candidate_chain_ids,assembly_id\n')
            rejected_csv.write(
                'comp_name,type,chain_ids,candidate_pdb_id,'
                'candidate_chain_ids,assembly_id,reason\n')

        for comp_name, candidates in processed_csv.items():
            if comp_name in post_processed:
                continue

            pdb_id, ab_chains, ag_chains = comp_name_to_pdb_and_chains(
                comp_name)

            ab_candidates_pdb_ids = get_pdb_ids(candidates, AB)
            ag_candidates_pdb_ids = get_pdb_ids(candidates, AG)

            filter_candidates_pack(comp_name, pdb_id, ab_candidates_pdb_ids,
                                   ab_chains, AB,
                                   filtered_csv, rejected_csv)

            filter_candidates_pack(comp_name, pdb_id, ag_candidates_pdb_ids,
                                   ag_chains, AG,
                                   filtered_csv, rejected_csv)

            post_processed_csv.write(str(comp_name) + '\n')
            post_processed_csv.flush()


if __name__ == '__main__':
    filter_for_unboundness(
        process_csv(pd.read_csv('unbound_data.csv').drop_duplicates()))
