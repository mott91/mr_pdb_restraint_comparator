import re
from pathlib import Path

def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def convert_atom_name(atom_name):
    # Remove any quotation marks
    atom_name = atom_name.strip('"')
    
    # Handle special cases for H5' and H5''
    if atom_name.startswith("H5'"):
        if atom_name.endswith('1'):
            return "H5'"
        elif atom_name.endswith('2'):
            return "H5''"
    
    # Handle other special cases (H2', H2'', etc.)
    if "'" in atom_name and atom_name[-1].isdigit():
        base_name = atom_name[:-1]
        if atom_name.endswith('1'):
            return base_name
        elif atom_name.endswith('2'):
            return base_name + "'"
    
    # For all other cases, return the original name
    return atom_name

def parse_mr_file(mr_file_path):
    mr_data = []
    header_block = []
    in_header = False

    try:
        with open(mr_file_path, 'r') as file:
            for line in file:
                line = line.strip()

                if line.startswith('!'):
                    if not in_header and line == '!':
                        in_header = True
                    header_block.append(line)
                    if in_header and line == '!':
                        mr_data.append(('header', header_block))
                        header_block = []
                        in_header = False
                elif line.startswith('assign'):
                    if header_block:
                        mr_data.append(('header', header_block))
                        header_block = []
                        in_header = False
                    parts = re.findall(r'\(resid\s+(\d+)\s+and\s+name\s+([\w\'\"]+)\s*\)', line)
                    parts = [(resid, convert_atom_name(atom_name.strip('"'))) for resid, atom_name in parts]
                    float_values = re.findall(r'\d+\.\d+', line)
                    mr_data.append(('data', (parts, float_values)))
                else:
                    if header_block:
                        mr_data.append(('header', header_block))
                        header_block = []
                        in_header = False

        if header_block:
            mr_data.append(('header', header_block))

    except FileNotFoundError:
        print(f"Error: The file {mr_file_path} was not found.")
    except Exception as e:
        print(f"An error occurred while parsing the MR file: {e}")

    return mr_data

def parse_pdb_file(pdb_file_path):
    pdb_data = {}

    try:
        with open(pdb_file_path, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    atom_index = int(line[6:11])
                    atom_name = convert_atom_name(line[12:16].strip())
                    resid = int(line[22:26])
                    pdb_data[(resid, atom_name)] = atom_index

    except FileNotFoundError:
        print(f"Error: The file {pdb_file_path} was not found.")
    except Exception as e:
        print(f"An error occurred while parsing the PDB file: {e}")

    return pdb_data

def link_mr_to_pdb(mr_data, pdb_data):
    linked_data = []
    warnings = {}

    for item in mr_data:
        if item[0] == 'header':
            linked_data.append(item)
        else:
            parts, float_values = item[1]
            indices = []
            for resid, atom_name in parts:
                key = (int(resid), atom_name)
                if key in pdb_data:
                    indices.append(str(pdb_data[key]))
                else:
                    print(f"Warning: Atom {atom_name} within resid {resid} not found in PDB data.")
                    indices.append("N/A")
                    if (resid, atom_name) in warnings:
                        warnings[(resid, atom_name)] += 1
                    else:
                        warnings[(resid, atom_name)] = 1
            atom_names = [atom_name for _, atom_name in parts]
            linked_data.append(('data', (indices, atom_names, float_values)))

    return linked_data, warnings

def format_output_line(indices, atom_names, float_values, include_atom_names=False):
    if include_atom_names:
        formatted_indices = [f"{index:>6}({atom_name})" for index, atom_name in zip(indices, atom_names)]
    else:
        formatted_indices = [f"{index:>6}" for index in indices]

    formatted_floats = [f"{float(value):<8.3f}" for value in float_values]
    separator = " " * 4

    return separator.join([" ".join(formatted_indices), " ".join(formatted_floats)]).rstrip()

def main(include_atom_names=False):
    mr_file_path = Path('/home/mot/14mer/8clr.mr')
    pdb_file_path = Path('/home/mot/14mer/8clr.pdb')
    output_file_path = mr_file_path.parent / 'nmr_restraints.txt'

    mr_data = parse_mr_file(mr_file_path)
    if not mr_data:
        print("Error: No MR data found.")
        return

    pdb_data = parse_pdb_file(pdb_file_path)
    if not pdb_data:
        print("Error: No PDB data found.")
        return

    linked_data, warnings = link_mr_to_pdb(mr_data, pdb_data)

    try:
        with open(output_file_path, 'w') as file:
            for item in linked_data:
                if item[0] == 'header':
                    file.writelines(f"{line.replace('!', '#')}\n" for line in item[1])
                else:
                    indices, atom_names, float_values = item[1]
                    output_line = format_output_line(indices, atom_names, float_values, include_atom_names)
                    file.write(f"{output_line}\n")

        print(f"Output has been written to {output_file_path}")

        # Print warning summary
        total_warnings = sum(warnings.values())
        unique_warnings = len(warnings)
        print(f"\nWarning Summary:")
        print(f"Total atom names not found: {total_warnings}")
        print(f"Unique atom name/residue combinations not found: {unique_warnings}")
        print("\nTop 5 most frequent warnings:")
        for (resid, atom_name), count in sorted(warnings.items(), key=lambda x: x[1], reverse=True)[:5]:
            print(f"  Residue {resid}, Atom {atom_name}: {count} times")

    except Exception as e:
        print(f"An error occurred while writing to the output file: {e}")

if __name__ == "__main__":
    main(include_atom_names=False)  # Set to True to include atom names in the output # useful for debugging
