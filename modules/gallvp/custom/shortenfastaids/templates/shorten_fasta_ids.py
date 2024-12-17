#!/usr/bin/env python3

import re
from importlib.metadata import version
from platform import python_version

from Bio import SeqIO

# The input fasta file path
fasta_file_path = "$fasta"
output_files_prefix = "$prefix"


def extract_fasta_ids_and_descriptions(fasta_file_path):
    fasta_file_obj = SeqIO.parse(fasta_file_path, "fasta")

    ids = []
    for record in fasta_file_obj:
        ids.append((record.id, record.description))
    return ids


def write_fasta_with_new_ids(fasta_file_path, id_mapping, file_prefix):
    old_fasta_file_obj = SeqIO.parse(fasta_file_path, "fasta")
    id_map = dict(id_mapping)

    replaced_records = []
    for record in old_fasta_file_obj:
        old_id = record.id

        new_id = id_map[old_id]
        record.id = new_id
        record.description = ""

        replaced_records.append(record)

    SeqIO.write(replaced_records, f"{file_prefix}.short.ids.fasta", "fasta")


def do_id_need_to_change(id_and_description, silent=False):
    id = id_and_description[0]
    description = id_and_description[1]
    if len(id) > 13:
        if not silent:
            print(f"{id} has length greater than 13")
        return True

    if not re.match(r"^[a-zA-Z0-9_]+\$", id):
        if not silent:
            print(f"{id} does not match '^[a-zA-Z0-9_]+\$'")
        return True

    if description != id and description != "":
        if not silent:
            print(f"{id} contains a comment: {description.replace(id, '')}")
        return True

    if not silent:
        print(f"{id} is acceptable")
    return False


def do_ids_need_to_change(ids_and_descriptions, silent=False):
    return any([do_id_need_to_change(id_and_description, silent) for id_and_description in ids_and_descriptions])


def extract_common_patterns(ids):
    pattern_counts = {}
    for id in ids:
        patterns = re.findall(r"[A-Za-z0_]{4,}", id)
        for pattern in set(patterns):
            pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1

    common_patterns = [pattern for pattern, count in pattern_counts.items() if count >= 2]

    if len(common_patterns) < 1:
        return {}

    return {pattern: pattern[:3] for pattern in common_patterns}


def shorten_ids(input_ids_and_descriptions, patterns_dict):
    shortened_ids = []

    for id_and_description in input_ids_and_descriptions:
        id = id_and_description[0]
        description = ""  # Treat description as absent as it will be removed by write_fasta_with_new_ids
        if not do_id_need_to_change((id, description), silent=True):
            shortened_ids.append(id)
            continue

        shortened_id = shorten_id_by_pattern_replacement(patterns_dict, id)

        if not do_id_need_to_change((shortened_id, description), silent=True):
            shortened_ids.append(shortened_id)
            continue

        shortened_id = f"Ctg{generate_hash(id)}"

        if not do_id_need_to_change((shortened_id, description), silent=True):
            shortened_ids.append(shortened_id)
            continue

        raise ValueError(f"Failed to shorten id: {id} ({shortened_id})")

    return shortened_ids


def shorten_id_by_pattern_replacement(patterns_dict, id):
    if patterns_dict == {}:
        return id

    shortened_id = id
    matches_for_id = match_substrings(patterns_dict.keys(), shortened_id)

    for pattern in matches_for_id:
        shortened_id = re.sub(
            rf"({re.escape(pattern)})",
            patterns_dict[pattern],
            shortened_id,
        )
    return shortened_id if shortened_id[len(shortened_id) - 1] != "_" else shortened_id[0 : (len(shortened_id) - 1)]


def match_substrings(substrings, target_string):
    pattern = "|".join(map(re.escape, substrings))
    matches = re.findall(pattern, target_string)
    return matches


def generate_hash(string):
    import hashlib

    hash_object = hashlib.sha1(string.encode())
    full_hash = hash_object.hexdigest()
    short_hash = full_hash[:10]
    return short_hash


def fail_if_new_ids_not_valid(ids):
    if len(ids) != len(set(ids)):
        raise ValueError("Th new IDs are not unique")


if __name__ == "__main__":
    input_ids_and_descriptions = extract_fasta_ids_and_descriptions(fasta_file_path)
    input_ids = [x[0] for x in input_ids_and_descriptions]

    # Write versions
    with open("versions.yml", "w") as f_versions:
        f_versions.write('"${task.process}":\\n')
        f_versions.write(f"    python: {python_version()}\\n")
        f_versions.write(f"    biopython: {version('biopython')}\\n")

    if not do_ids_need_to_change(input_ids_and_descriptions):
        print("IDs have acceptable length and character. No change required.")
        with open(f"{output_files_prefix}.short.ids.tsv", "w") as f:
            f.write("IDs have acceptable length and character. No change required.")
        exit(0)

    new_ids = shorten_ids(input_ids_and_descriptions, extract_common_patterns(input_ids))
    fail_if_new_ids_not_valid(new_ids)

    with open(f"{output_files_prefix}.short.ids.tsv", "w") as f:
        for input_id, new_id in zip(input_ids, new_ids):
            f.write(f"{input_id}\\t{new_id}\\n")

    write_fasta_with_new_ids(fasta_file_path, zip(input_ids, new_ids), output_files_prefix)
