import os
import argparse

def count_atoms_in_pdb(pdb_file):
    """Count number of atoms in a PDB file."""
    count = 0
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                count += 1
    return count

def calculate_density(pdb_file, out_file):
    """Calculate density using atom count from PDB and last line of .out file."""
    # Get atom count from pdb
    atom_count = count_atoms_in_pdb(pdb_file)

    if not os.path.exists(out_file):
        raise FileNotFoundError(f"File not found: {out_file}")

    with open(out_file, 'r') as file:
        lines = file.readlines()
        if not lines:
            raise ValueError(f"No data in {out_file}")

        # Take the last line and parse
        last_line = lines[-1].strip()
        try:
            second_column_value = float(last_line.split()[1])
        except (IndexError, ValueError):
            raise ValueError(f"Could not parse numeric value from last line of {out_file}")

        # Formula: result = atom_count * 1000 / second_column_value
        density = atom_count * 1000 / second_column_value # res/nm^3
        return density

def main():
    parser = argparse.ArgumentParser(description="Calculate density from PDB and OUT file")
    parser.add_argument("pdb_file", help="Path to input PDB file")
    parser.add_argument("out_file", help="Path to corresponding .out file")

    args = parser.parse_args()

    density = calculate_density(args.pdb_file, args.out_file)
    print(f"Density = {density:.4f} res/nm^3")

if __name__ == "__main__":
    main()
