import math
#we open the pdb file
pdb_file = "/Users/black/Downloads/7u6r.pdb"
#we define a function, it will call pdbfile, and the user will input two residues, int the form eg 1, A
def calculate_distance(pdb_file, residue1, residue2):
    #we create a tuple where the values the user gave are given as the res number and the chain letter
    #for each residue
    resnum1, chain1 = residue1
    resnum2, chain2 = residue2

    #we create an empty variable to store the coodinates and to define the range
    coords1 = None
    coords2 = None

    first_resnumA = None
    last_resnumA = None

    first_resnumB = None
    last_resnumB = None

    #we create a with loop where it opens the pdb in read mode, for a given file, f. for each line in f:
    with open(pdb_file, 'r') as f:
        for line in f:
            chain = line[21].strip()
            resnum = line[22:26].strip()
            atom_name = line[12:16].strip()

            #how many resiudes are in chain A?
            if atom_name == "CA" and chain == "A":
                if first_resnumA is None:
                    first_resnumA = resnum
                last_resnumA = resnum


            if atom_name == "CA" and chain == "B":
                if first_resnumB is None:
                    first_resnumB = resnum
                last_resnumB = resnum

            if line.startswith("ATOM") and atom_name == "CA" and chain == chain1 and resnum == str(resnum1):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords1 =(x, y, z)

            if line.startswith("ATOM") and atom_name == "CA" and chain == chain2 and resnum == str(resnum2):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords2 = (x, y, z)

    print(f"Chain A resiudes range from {first_resnumA} to {last_resnumA}")
    print(f"Chain B residues range from {first_resnumB} to {last_resnumB}")

    if coords1 is None or coords2 is None:
        return "Error: One of both residues not found."

    distance = math.sqrt(
        (coords2[0] - coords1[0]) ** 2 + (coords2[1] - coords1[1]) ** 2 + (coords2[2] - coords1[2]) ** 2)
    return distance

result = calculate_distance(pdb_file, (36, "A"), (38    , "A"))
print(f"Distance: {result} Angstroms")

