import argparse

# Residue classification dictionary.
type_dict = {
    'ALA': 'N', 'ARG': 'B', 'ASN': 'P', 'ASP': 'A', 'CYS': 'A', 'GLU': 'A',
    'GLN': 'P', 'GLY': 'G', 'HIS': 'B', 'ILE': 'N', 'LEU': 'N', 'LYS': 'B',
    'MET': 'N', 'PHE': 'N', 'PRO': 'N', 'SER': 'P', 'THR': 'P', 'TRP': 'N',
    'TYR': 'A', 'VAL': 'N', 'NTR': 'B', 'CTR': 'A',
    'ACE': 'C', 'NME': 'C'
}

# One-letter to three-letter residue conversion.
one_to_three = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

def read_fasta_concat(fasta_file):
    sequence = ""
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                sequence += line.upper()
    return sequence


def generate_charge_file(fasta_file, output_file="charges.txt"):

    seq = read_fasta_concat(fasta_file)

    with open(output_file, "w") as out:
        for i, aa in enumerate(seq):
            if aa not in one_to_three:
                raise ValueError(
                    f"Unsupported residue at position {i + 1}: '{aa}'"
                )

            res3 = one_to_three[aa]
            tipo = type_dict[res3]

            if tipo == 'A':
                charge = -1.0
            elif tipo == 'B':
                charge = 1.0
            else:
                charge = 0.0

            out.write(f"{i} {charge}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a charge file from a FASTA sequence."
    )
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument(
        "-o",
        "--output",
        default="charges.txt",
        help="Output charge file",
    )
    args = parser.parse_args()

    generate_charge_file(args.fasta_file, args.output)

if __name__ == "__main__":
    main()
