#!/usr/bin/env python3

from platform import python_version

ids_tsv = "$ids_tsv"
input_gff3 = "$gff3"
output_prefix = "$prefix"


def create_name_mapping_from_tsv(file_path):
    dictionary = {}

    with open(file_path) as tsv_file:
        for line in tsv_file:
            columns = line.strip().split("\\t")
            if len(columns) != 2:
                raise ValueError(f"{file_path} should be a two column TSV file")

            orig_id, new_id = columns[0], columns[1]
            dictionary[new_id] = orig_id

    return dictionary


def restore_gff3_ids(new_to_orig_ids, file_path, output_file_name):
    # Write versions
    with open("versions.yml", "w") as f_versions:
        f_versions.write('"${task.process}":\\n')
        f_versions.write(f"    python: {python_version()}\\n")

    with open(file_path) as input_gff3_file:
        input_lines = input_gff3_file.readlines()

    with open(output_file_name, "w") as output_gff_file:
        for line in input_lines:
            if line.startswith("##"):
                output_gff_file.write(line)
                continue

            new_id = line.split("\\t")[0]
            orig_id = new_to_orig_ids[new_id]
            output_gff_file.write("\\t".join([orig_id] + line.split("\\t")[1:]))


if __name__ == "__main__":
    new_to_orig_ids = create_name_mapping_from_tsv(ids_tsv)
    restore_gff3_ids(new_to_orig_ids, input_gff3, f"{output_prefix}.restored.ids.gff3")
