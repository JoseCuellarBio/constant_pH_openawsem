import numpy as np
import pandas as pd
import math
import random
import os


class Protein:
    def __init__(self, source, coords=None, list_charged_residues=None, pH=7.0):

        self.type_dict = {
            'ALA': 'N', 'ARG': 'B', 'ASN': 'P', 'ASP': 'A', 'CYS': 'A', 'GLU': 'A',
            'GLN': 'P', 'GLY': 'G', 'HIS': 'B', 'ILE': 'N', 'LEU': 'N', 'LYS': 'B',
            'MET': 'N', 'PHE': 'N', 'PRO': 'N', 'SER': 'P', 'THR': 'P', 'TRP': 'N',
            'TYR': 'A', 'VAL': 'N', 'NTR': 'B', 'CTR': 'A',
            'ACE': 'C', 'NME': 'C'
        }

        self.list_charged_residues = list_charged_residues or []

        # Detect input type.
        if isinstance(source, str):
            self.data = self._parse_pdb(source)

        elif isinstance(source, list) and coords is not None:
            # OpenAWSEM input.
            self.data = self._parse_openawsem(source, coords)

        else:
            # DataFrame input.
            self.data = self._parse_dataframe(source)

        if not self.list_charged_residues:
            self._generate_provisional_charges()

        self.mc = MonteCarloResidue(self)
        self.neighborhood = Neighborhood(self)
        self.protonation_mc = ProtonationMC(self, pH=pH)

    # -----------------------------------

    def _generate_provisional_charges(self):
        charges = []

        for resid, info in self.data.items():

            if info["type"] == "A":
                info["charge"] = -1.0
                charges.append((resid, -1.0))

            elif info["type"] == "B":
                info["charge"] = 1.0
                charges.append((resid, 1.0))

            else:
                info["charge"] = 0.0

        self.list_charged_residues = charges

    # -----------------------------------

    def _parse_openawsem(self, target_atoms_info, coords):

        data = {}
        charge_map = dict(self.list_charged_residues)

        for (atom_index, atom_name, resid, resname), pos in zip(target_atoms_info, coords):

            x, y, z = pos

            if resid not in data:

                data[resid] = {
                    "resname": resname,
                    "type": self.type_dict.get(resname, "X"),
                    "charge": charge_map.get(resid, 0.0),
                    "atoms": {}
                }

            data[resid]["atoms"][atom_name] = {
                "coords": [float(x), float(y), float(z)],
                "index": atom_index
            }

        return data
    # -----------------------------------

    def _parse_dataframe(self, df):

        data = {}
        charge_map = dict(self.list_charged_residues)

        for _, row in df.iterrows():

            resid = int(row["residue_number"])
            resname = row["residue_name"]
            atom = row["atom"]
            x, y, z = row["posicion_xyz"]

            if resid not in data:

                data[resid] = {
                    "resname": resname,
                    "type": self.type_dict.get(resname, "X"),
                    "charge": charge_map.get(resid, 0.0),
                    "atoms": {}
                }

            data[resid]["atoms"][atom] = {
                "coords": [float(x), float(y), float(z)]
            }

        return data

class MonteCarloResidue:
    def __init__(self, protein):
        self.protein = protein
        self.current_resid = None

    def _get_charged_resids(self):
        return [resid for resid, _ in self.protein.list_charged_residues]

    def choose_residue(self):
        charged = self._get_charged_resids()
        self.current_resid = random.choice(charged) if charged else None
        return self.current_resid

    def get_current_residue_info(self):
        if self.current_resid is None:
            return None
        return {self.current_resid: self.protein.data[self.current_resid]}

class Neighborhood:
    def __init__(self, protein, cutoff=6.0):
        self.protein = protein
        self.cutoff = cutoff

    def _get_atoms_neighborhood(self, center_resid):
        coords_list, resids = [], []

        for resid, info in self.protein.data.items():
            if resid == center_resid:
                continue

            atom = next((a for a in ["CB", "CA", "O"] if a in info["atoms"]), None)
            if atom:
                coords_list.append(info["atoms"][atom]["coords"])
                resids.append(resid)

        return np.array(coords_list), resids

    def compute_distances(self):
        center_resid = self.protein.mc.current_resid
        if center_resid is None:
            return None

        center_atoms = self.protein.data[center_resid]["atoms"]
        atom = next((a for a in ["CB", "CA", "O"] if a in center_atoms), None)
        center_coords = np.array(center_atoms[atom]["coords"])

        neighbor_coords, resids = self._get_atoms_neighborhood(center_resid)
        distances = np.linalg.norm(neighbor_coords - center_coords, axis=1)
        return resids, distances

    def dataframe(self):
        result = self.compute_distances()
        if result is None:
            return None

        resids, distances = result
        return pd.DataFrame({
            "resid": resids,
            "resname": [self.protein.data[r]["resname"] for r in resids],
            "type": [self.protein.data[r]["type"] for r in resids],
            "charge": [self.protein.data[r]["charge"] for r in resids],
            "distance": distances
        })

def calculate_new_charge(resid, acid_basic, charge):
    if acid_basic == -1:
        new_charge = 0.0 if charge == -1 else -1.0
    elif acid_basic == 1:
        new_charge = 1.0 if charge == 0 else 0.0
    else:
        new_charge = charge
    return new_charge, new_charge - charge

def calculate_delta_term_elec(df, residue_mc, new_charge):
    K_elec = 2.43232 / 10
    L = 1.0
    resid_mc = list(residue_mc.keys())[0]
    old_charge = residue_mc[resid_mc]["charge"]

    def calc(qc):
        e = 0
        for _, row in df.iterrows():
            qn, r = row["charge"], row["distance"]
            if r == 0 or qn in [None, 0]:
                continue
            sign = 1 if qn > 0 else -1
            e += qc * ((sign / r) * math.exp(-r / L))
        return e * K_elec

    return calc(new_charge) - calc(old_charge)

def calculate_delta_term_pH(resname, delta_q, pH, T=300):
    pKa = {'ASP':4.0,'GLU':4.5,'LYS':10.6,'ARG':12.0,'HIS':6.4,'CYS':8.3,'TYR':11.0,'NTR':7.5,'CTR':3.5}
    kb = 0.001987
    return delta_q * (pH - pKa[resname]) * kb * T * np.log(10)

def accept_or_reject(resid, protein, charged_residues,
                     dE_pH, dE_elec, dE_polar, new_charge):

    kb, T = 0.001987, 300
    dE = dE_pH + dE_elec + dE_polar
    # print("deltaE", dE)

    new_list = charged_residues.copy()

    if dE < 0 or random.random() < math.exp(-dE/(kb*T)):

        # Update the charged-residue list.
        new_list = [(r, new_charge if r == resid else q) for r, q in charged_residues]

        # Get the particle index (CB or CA).
        atom_dict = protein.data[resid]["atoms"]
        atom_name = list(atom_dict.keys())[0]
        particle_index = atom_dict[atom_name]["index"]

        return new_list, (particle_index, new_charge)

    else:
        return charged_residues, None

class ProtonationMC:

    def __init__(self, protein, pH=7.0, T=300):
        self.protein = protein
        self.pH = pH
        self.T = T


    def attempt_charge_flip(self, charged_residues):

        df = self.protein.neighborhood.dataframe()
        residue_mc = self.protein.mc.get_current_residue_info()

        if df is None or residue_mc is None:
            return charged_residues, None

        resid = list(residue_mc.keys())[0]
        info = residue_mc[resid]

        old_charge = info["charge"]

        # print("DEBUG residue type:", resid, info["resname"], "->", info["type"])

        acid_base = -1 if info["type"] == "A" else 1 if info["type"] == "B" else 0

        if acid_base == 0:
            return charged_residues, None

        new_charge, delta_q = calculate_new_charge(resid, acid_base, old_charge)

        dE_elec = calculate_delta_term_elec(df, residue_mc, new_charge)
        dE_polar = 0.0
        dE_pH = calculate_delta_term_pH(info["resname"], delta_q, self.pH, self.T)

        new_list, particle_info = accept_or_reject(
            resid,
            self.protein,
            charged_residues,
            dE_pH,
            dE_elec,
            dE_polar,
            new_charge
        )

        accepted = new_list != charged_residues

        if accepted:
            self.protein.data[resid]["charge"] = new_charge

        return new_list, particle_info
        
#############################################
##### MORE FUNCTIONS ########################
#############################################

def process_charged_residue_file(filename):
    with open(filename, 'r') as file:
        residues = []
        for line in file:
            parts = line.split()
            residue = int(parts[0])
            charge = float(parts[1])
            residues.append((residue, charge))

    # Keep residues with nonzero charge.
    charged_residues = [(residue, charge) for residue, charge in residues if charge != 0.0]
    
    return charged_residues

def get_target_atom_indices_and_info(oa):  
    """Return target indices and full CB/CA atom metadata for efficient filtering."""  
    target_indices = []  
    target_atoms_info = []  
      
    one_to_three = {  
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP",  
        "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",  
        "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",  
        "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",  
        "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"  
    }  
      
    for residue in oa.pdb.topology.residues():  
        cb_atom = None  
        ca_atom = None  
          
        for atom in residue.atoms():  
            if atom.name == 'CB':  
                cb_atom = atom  
            elif atom.name == 'CA':  
                ca_atom = atom  
          
        target_atom = cb_atom if cb_atom is not None else ca_atom  
          
        if target_atom is not None:  
            # Index for fast position lookup.
            target_indices.append(target_atom.index)  
              
            # Full residue information.
            if residue.index < len(oa.seq):  
                real_resname_one = oa.seq[residue.index]  
                real_resname_three = one_to_three.get(real_resname_one, "UNK")  
            else:  
                real_resname_three = "UNK"  
              
            target_atoms_info.append(  
                (target_atom.index, target_atom.name, residue.index, real_resname_three)  
            )  
      
    return target_indices, target_atoms_info
