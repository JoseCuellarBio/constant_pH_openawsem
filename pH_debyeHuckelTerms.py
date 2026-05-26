try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ModuleNotFoundError:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
import numpy as np
    
def debye_huckel_term_constant_ph(oa, k_dh=5*4.184, forceGroup=30, screening_length=1.0, chargeFile=None):

    print("Debye Huckel term is ON")

    k_dh *= oa.k_awsem * 0.1
    k_screening = 1.0

    #dh = CustomNonbondedForce(
    #    f"{k_dh}*charge1*charge2/r*exp(-{k_screening}*r/{screening_length})"
    #    "*step(abs(res1-res2)-2+isChainEdge1*isChainEdge2)"
    #)
    
    dh = CustomNonbondedForce(
        f"{k_dh}*charge1*charge2/r*exp(-{k_screening}*r/{screening_length})"
        "*step(abs(res1-res2)+isChainEdge1*isChainEdge2)"
    )
    
    dh.addPerParticleParameter("charge")
    dh.addPerParticleParameter("res")
    dh.addPerParticleParameter("isChainEdge")

    if oa.periodic:
        dh.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    else:
        dh.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)

    dh.setCutoffDistance(3.0)

    # Read charge file.

    if chargeFile is None:
        raise ValueError("chargeFile must be provided")

    chargeInfo = np.loadtxt(chargeFile)
    charge_dict = {int(i): float(q) for i, q in chargeInfo}

    # Use CB, or CA for GLY.

    cb_fixed = [x if x > 0 else y for x, y in zip(oa.cb, oa.ca)]

    charged_atoms = []

    # Assign charges.

    for i in range(oa.natoms):

        res_idx = oa.resi[i]
        charge = 0.0

        if res_idx in charge_dict and i in cb_fixed:
            charge = charge_dict[res_idx]

        if i in oa.chain_starts and i in oa.n:
           charge = 1.0
        
        if i in oa.chain_ends and i in oa.o:
           charge = -1.0
        
        is_chain_edge = 1 if (i in oa.chain_ends or i in oa.chain_starts) else 0
        dh.addParticle([charge, res_idx, is_chain_edge])

        if charge != 0.0:
            charged_atoms.append(i)

    # Restrict interactions to charged atoms.

    if len(charged_atoms) > 0:
        dh.addInteractionGroup(charged_atoms, charged_atoms)

    dh.setForceGroup(forceGroup)

    return dh
