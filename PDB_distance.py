import math

# Path to PDB file
pdb_file = "/Users/black/Downloads/7u6r.pdb"


def calculate_distance(pdb_file, residue1, residue2):
    """
    Calculate the distance between CA atoms of two residues in a PDB file.

    Args:
        pdb_file: Path to PDB file
        residue1: Tuple of (residue_number, chain) e.g., (36, "A")
        residue2: Tuple of (residue_number, chain) e.g., (37, "A")

    Returns:
        Distance in Angstroms, or error message if residues not found
    """

    # Extract residue information from input tuples
    resnum1, chain1 = residue1
    resnum2, chain2 = residue2

    # Initialize variables to store coordinates
    coords1 = None
    coords2 = None

    # Track residue ranges for each chain
    firstresA = None
    lastresA = None
    firstresB = None
    lastresB = None

    # Open and parse PDB file
    with open(pdb_file, "r") as f:
        for line in f:
            # Extract data from PDB format (fixed column positions)
            chain = line[21].strip()
            resnum = line[22:26].strip()
            atom_name = line[12:16].strip()

            # Track residue range for chain A
            if atom_name == "CA" and chain == "A":
                if firstresA is None:
                    firstresA = resnum
                lastresA = resnum

            # Track residue range for chain B
            if atom_name == "CA" and chain == "B":
                if firstresB is None:
                    firstresB = resnum
                lastresB = resnum

            # Extract coordinates for first residue
            if (resnum == str(resnum1) and chain == chain1 and
                    atom_name == "CA" and line.startswith("ATOM")):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords1 = (x, y, z)

            # Extract coordinates for second residue
            if (resnum == str(resnum2) and chain == chain2 and
                    atom_name == "CA" and line.startswith("ATOM")):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords2 = (x, y, z)

    # Display available residue ranges
    print(f"Chain A residues range from {firstresA} to {lastresA}")
    print(f"Chain B residues range from {firstresB} to {lastresB}")

    # Check if both residues were found
    if coords1 is None or coords2 is None:
        return "Error: One or both residues not found."

    # Calculate Euclidean distance
    distance = math.sqrt(
        (coords2[0] - coords1[0]) ** 2 +
        (coords2[1] - coords1[1]) ** 2 +
        (coords2[2] - coords1[2]) ** 2
    )

    return distance


# Example usage
result = calculate_distance(pdb_file, (36, "A"), (37, "A"))
print(f"Distance: {result} Angstroms")