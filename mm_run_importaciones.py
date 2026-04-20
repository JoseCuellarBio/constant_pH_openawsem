#!/usr/bin/env python3
import os
import sys
import random
import time
from random import seed, randint
import argparse
import platform
from datetime import datetime
from time import sleep
import fileinput
import importlib.util

from openawsem import *
from openawsem.helperFunctions.myFunctions import *
###change###
from openmm.app import PDBFile, PDBReporter
from openmm import Vec3
import json
from Montecarlo import *

do = os.system
cd = os.chdir


def run(args):
    simulation_platform = args.platform
    platform = Platform.getPlatformByName(simulation_platform)
    if simulation_platform == "CPU":
        if args.thread != -1:
            platform.setPropertyDefaultValue("Threads", str(args.thread))
        print(f"{simulation_platform}: {platform.getPropertyDefaultValue('Threads')} threads")

    # if mm_run.py is not at the same location of your setup folder.
    setupFolderPath = os.path.dirname(args.protein)
    setupFolderPath = "." if setupFolderPath == "" else setupFolderPath
    proteinName = pdb_id = os.path.basename(args.protein)


    pwd = os.getcwd()
    toPath = os.path.abspath(args.to)
    checkPointPath = None if args.fromCheckPoint is None else os.path.abspath(args.fromCheckPoint)
    forceSetupFile = None if args.forces is None else os.path.abspath(args.forces)
    parametersLocation = "." if args.parameters is None else os.path.abspath(args.parameters)
    os.chdir(setupFolderPath)

    # chain=args.chain.upper()
    chain=args.chain
    pdb = f"{pdb_id}.pdb"

    if chain == "-1":
        chain = getAllChains("crystal_structure.pdb")
        print("Chains to simulate: ", chain)

    if args.to != "./":
        # os.system(f"mkdir -p {args.to}")
        os.makedirs(toPath, exist_ok=True)
        os.system(f"cp {forceSetupFile} {toPath}/forces_setup.py")
        os.system(f"cp crystal_structure.fasta {toPath}/")
        os.system(f"cp crystal_structure.pdb {toPath}/")
        # os.system(f"cp {pdb} {args.to}/{pdb}")
        # pdb = os.path.join(args.to, pdb)

    if args.fromOpenMMPDB:
        input_pdb_filename = proteinName
        seq=read_fasta("crystal_structure.fasta")
        print(f"Using Seq:\n{seq}")
    else:
        suffix = '-openmmawsem.pdb'
        if pdb_id[-len(suffix):] == suffix:
            input_pdb_filename = pdb_id
        else:
            input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"
        seq=None

    if args.fasta == "":
        seq = None
    else:
        seq = seq=read_fasta(args.fasta)
        print(f"Using Seq:\n{seq}")
    # start simulation
    collision_rate = 5.0 / picoseconds
    checkpoint_file = "checkpnt.chk"
    checkpoint_reporter_frequency = 10000
    snapShotCount = 400
    stepsPerT = int(args.steps/snapShotCount)
    Tstart = args.tempStart
    Tend = args.tempEnd
    if args.reportFrequency == -1:
        if stepsPerT == 0:
            reporter_frequency = 4000
        else:
            reporter_frequency = stepsPerT
    else:
        reporter_frequency = args.reportFrequency

    print(f"using force setup file from {forceSetupFile}")
    spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
    forces = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(forces)


    oa = OpenMMAWSEMSystem(input_pdb_filename, k_awsem=1.0, chains=chain, xml_filename=openawsem.xml, seqFromPdb=seq, includeLigands=args.includeLigands, periodic=True)  # k_awsem is an overall scaling factor that will affect the relevant temperature scales
    box_size = 27.1 * nanometer

    oa.system.setDefaultPeriodicBoxVectors(
        Vec3(box_size, 0, 0),
        Vec3(0, box_size, 0),
        Vec3(0, 0, box_size)
    )
    
    myForces = forces.set_up_forces(oa, submode=args.subMode, contactParameterLocation=parametersLocation)
    oa.addForcesWithDefaultForceGroup(myForces)
    
    print("\n===== FULL PBC CHECK =====")
    for i, force in enumerate(oa.system.getForces()):
        print(i, type(force))
        try:
            print("   usesPBC:", force.usesPeriodicBoundaryConditions())
        except:
            print("   usesPBC: method not available")

    # Si es NonbondedForce mostramos método
    from openmm import NonbondedForce
    if isinstance(force, NonbondedForce):
        print("   NonbondedMethod:", force.getNonbondedMethod())
        print("   Cutoff:", force.getCutoffDistance())

    print("==========================\n")

    if args.fromCheckPoint:
        integrator = LangevinIntegrator(Tstart*kelvin, 1/picosecond, args.timeStep*femtoseconds)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        simulation.loadCheckpoint(checkPointPath)
    else:
        # output the native and the structure after minimization
        integrator = CustomIntegrator(0.001)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
        simulation.reporters.append(PDBReporter(os.path.join(toPath, "native.pdb"), 1))
        simulation.reporters.append(DCDReporter(os.path.join(toPath, "movie.dcd"), 1))
        simulation.step(int(1))
        simulation.minimizeEnergy()  # first, minimize the energy to a local minimum
        simulation.step(int(1))

        integrator = LangevinIntegrator(Tstart*kelvin, 1/picosecond, args.timeStep*femtoseconds)
        simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)
        simulation.context.setPositions(oa.pdb.positions)  # set the initial positions of the atoms
        simulation.context.setVelocitiesToTemperature(Tstart*kelvin)  # set the initial velocities of the atoms
        simulation.minimizeEnergy()  # first, minimize the energy to a local minimum
        
    print("reporter_frequency", reporter_frequency)
    simulation.reporters.append(StateDataReporter(stdout, reporter_frequency, step=True, potentialEnergy=True, temperature=True, volume=True))  # output energy and temperature during simulation
    simulation.reporters.append(StateDataReporter(os.path.join(toPath, "output.log"), reporter_frequency, step=True, potentialEnergy=True, temperature=True, volume=True)) # output energy and temperature to a file
    simulation.reporters.append(PDBReporter(os.path.join(toPath, "movie.pdb"), reportInterval=reporter_frequency))  # output PDBs of simulated structures
    simulation.reporters.append(DCDReporter(os.path.join(toPath, "movie.dcd"), reportInterval=reporter_frequency, append=True))
    print("Simulation Starts")
    start_time = time.time()

    last_force = myForces[-1]

    if args.simulation_mode == 0:
        
        interrupt_frequency = args.interruptFrequency

        # Obtener ambas listas  
        indices, target_atoms_info = get_target_atom_indices_and_info(oa)  
         
        charged_residues = procesador_de_archivo_con_residuos_cargados('charge.txt')  # lee el archivo y los pasa a una lista de tuplas (Residuo,carga)

        seq_oa = oa.seq
        print("seq_oa", seq_oa)

        Hawsem_state = []

        # limpiar/crear archivo de estado
        with open('Hawsem.state', 'w') as f:
            pass

        # pH fijo para toda la simulación
        pH = 2

        total_steps = 0

        while total_steps < int(args.steps):

            # print("interrupt frequency", interrupt_frequency)

            remaining_steps = min(interrupt_frequency, int(args.steps) - total_steps)

            simulation.step(remaining_steps)

            total_steps += remaining_steps

            # Obtener posiciones actuales  
            state = simulation.context.getState(getPositions=True)  
            positions = state.getPositions()  
            coords = [(positions[i].x, positions[i].y, positions[i].z) for i in indices]

            # crear objeto proteína
            prot = Proteina(target_atoms_info, coords=coords, list_charged_residues=charged_residues, pH=pH)

            # elegir residuo para MC
            residue_mc = prot.mc.choose_residue()
            # print(residue_mc)

            # intento de cambio de protonación
            charged_residues, new_parameters = prot.protonation_mc.attempt_charge_flip(charged_residues)

            prot.list_charged_residues = charged_residues

            # actualizar parámetros de Debye-Huckel si hubo cambio
            if new_parameters is not None:

                particle_index, new_charge = new_parameters

                last_force.setParticleParameters(
                    particle_index,
                    [new_charge]
                )

                last_force.updateParametersInContext(simulation.context)

            # guardar estado de cargas
            if total_steps % reporter_frequency == 0:

                Hawsem_state_line = json.dumps(charged_residues) + str(pH)

                with open('Hawsem.state', 'a') as f:
                    f.write(Hawsem_state_line + "\n")

###Change####

    elif args.simulation_mode == 1:
        deltaT = (Tend - Tstart) / snapShotCount
        for i in range(snapShotCount):
            integrator.setTemperature((Tstart + deltaT*i)*kelvin)
            simulation.step(stepsPerT)

    time_taken = time.time() - start_time  # time_taken is in seconds
    hours, rest = divmod(time_taken,3600)
    minutes, seconds = divmod(rest, 60)
    print(f"---{hours} hours {minutes} minutes {seconds} seconds ---")

    timeFile = os.path.join(toPath, "time.dat")
    with open(timeFile, "w") as out:
        out.write(str(time_taken)+"\n")

    # accompany with analysis run
    simulation = None
    time.sleep(10)
    os.chdir(pwd)
    print(os.getcwd())
    if args.fasta == "":
        analysis_fasta = ""
    else:
        analysis_fasta = f"--fasta {args.fasta}"
    if args.includeLigands:
        additional_cmd = "--includeLigands"
    else:
        additional_cmd = ""
    os.system(f"{sys.executable} mm_analyze.py {args.protein} -t {os.path.join(toPath, 'movie.dcd')} --subMode {args.subMode} -f {args.forces} {analysis_fasta} {additional_cmd} -c {chain}")

def main():
    # from run_parameter import *
    parser = argparse.ArgumentParser(
        description="This is a python3 script to\
        automatic copy the template file, \
        run simulations")

    parser.add_argument("protein", help="The name of the protein")
    parser.add_argument("--name", default="simulation", help="Name of the simulation")
    parser.add_argument("--to", default="./", help="location of movie file")
    parser.add_argument("-c", "--chain", type=str, default="-1")
    parser.add_argument("-t", "--thread", type=int, default=-1, help="default is using all that is available")
    parser.add_argument("-p", "--platform", type=str, default="OpenCL")
    parser.add_argument("-s", "--steps", type=float, default=2e4, help="step size, default 1e5")
    parser.add_argument("--tempStart", type=float, default=800, help="Starting temperature")
    parser.add_argument("--tempEnd", type=float, default=200, help="Ending temperature")
    parser.add_argument("--fromCheckPoint", type=str, default=None, help="The checkpoint file you want to start from")
    parser.add_argument("-m", "--simulation_mode", type=int, default=1,
                    help="default 1,\
                            0: constant temperature,\
                            1: temperature annealing")
    parser.add_argument("--subMode", type=int, default=-1)
    parser.add_argument("-f", "--forces", default="forces_setup.py")
    parser.add_argument("--parameters", default=None)
    parser.add_argument("-r", "--reportFrequency", type=int, default=-1, help="default value step/400")
    parser.add_argument("--fromOpenMMPDB", action="store_true", default=False)
    parser.add_argument("--fasta", type=str, default="crystal_structure.fasta")
    parser.add_argument("--timeStep", type=float, default=2)
    parser.add_argument("--includeLigands", action="store_true", default=False)
    parser.add_argument("--interruptFrequency", type=int, default=1000, help="Frequency of interruptions during simulation")
    args = parser.parse_args()


    with open('commandline_args.txt', 'a') as f:
        f.write(' '.join(sys.argv))
        f.write('\n')
    print(' '.join(sys.argv))

    run(args)

if __name__=="__main__":
    main()
