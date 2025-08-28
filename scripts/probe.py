import os
import math
import argparse
def parse_pdb(file_path):
    """Extract atom coordinates from a PDB file."""
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append((x, y, z))
    return atoms

def calculate_center_of_mass(atoms):
    """Calculate center of mass from a list of atom coordinates."""
    if not atoms:
        return (0.0, 0.0, 0.0)
    x_sum = sum(x for x, _, _ in atoms)
    y_sum = sum(y for _, y, _ in atoms)
    z_sum = sum(z for _, _, z in atoms)
    num_atoms = len(atoms)
    return (x_sum / num_atoms, y_sum / num_atoms, z_sum / num_atoms)

def group_atoms_into_chains(atoms, chain_size):
    """Group atoms into chains of a given size."""
    return {
        f"Chain_{i // chain_size + 1}": atoms[i:i + chain_size]
        for i in range(0, len(atoms), chain_size)
    }

def calculate_average_distance_to_com(atoms, center_of_mass):
    """Calculate average distance of atoms in a chain to the cluster COM."""
    if not atoms:
        return 0.0
    return sum(
        math.dist(atom, center_of_mass)
        for atom in atoms
    ) / len(atoms)

def calculate_distance(point1, point2):
    """Euclidean distance between two 3D points."""
    return math.dist(point1, point2)


# ------------------------- MAIN ANALYSIS -------------------------


def main():
    parser = argparse.ArgumentParser(description="SPACEBALL_DENSITY – Surface chain analysis")
    parser.add_argument("pdb_file", help="Path to input PDB file")
    parser.add_argument("--chain_size", type=int, required=True, help="Number of atoms per chain")
    parser.add_argument("--surface_fractions", type=float, nargs="+", default=[0.5],
                        help="Surface fractions to analyze (e.g. 0.2 0.4 0.5)")

    args = parser.parse_args()

    if not os.path.exists(args.pdb_file):
        print(f"File not found: {args.pdb_file}")
        return

    # Parse structure
    atoms = parse_pdb(args.pdb_file)
    chains = group_atoms_into_chains(atoms, args.chain_size)
    total_chains = len(chains)
    cluster_com = calculate_center_of_mass(atoms)

    # Compute average distance of each chain to the cluster COM
    chain_distances = {
        chain_id: calculate_average_distance_to_com(chain_atoms, cluster_com)
        for chain_id, chain_atoms in chains.items()
    }

    # Sort chains by distance to COM (descending → surface chains first)
    sorted_chains = sorted(chain_distances.items(), key=lambda x: x[1], reverse=True)

    for frac in args.surface_fractions:
        N = max(1, int(total_chains * frac))
        top_surface_chains = sorted_chains[:N]

        print(f"\nTotal chains: {total_chains}, Surface fraction: {frac:.2f} → {N} chains")

        # Calculate COMs of the selected surface chains
        surface_coms = {
            chain_id: calculate_center_of_mass(chains[chain_id])
            for chain_id, _ in top_surface_chains
        }

        # Calculate minimum distance from each chain to another
        min_distances = []
        chain_ids = list(surface_coms.keys())

        for i, chain_id in enumerate(chain_ids):
            com1 = surface_coms[chain_id]
            min_distance = min(
                calculate_distance(com1, surface_coms[other_id])
                for j, other_id in enumerate(chain_ids) if i != j
            )
            min_distances.append(min_distance)
            #print(f"  {chain_id}: Min Distance = {min_distance:.2f}")

        # Output result
        if min_distances:
            avg_min_dist = sum(min_distances) / len(min_distances)
            print(f"  Probe_radius = {avg_min_dist:.2f} Å")
        else:
            print("  No valid distances computed.")


if __name__ == "__main__":
    main()
