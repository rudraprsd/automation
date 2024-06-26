import sys
import numpy as np
import os
import shutil
import subprocess
import time
from pathlib import Path
import yaml
import re
import pkg_resources

from collections import OrderedDict

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.inputs import Kpoints, Incar, Potcar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter
import matplotlib.pyplot as plt

def submit_job(script_path):
    """Submit the job to the queue."""
    try:
        result = subprocess.run(['qsub', script_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        if result.returncode != 0:
            raise RuntimeError(f"Job submission failed: {result.stderr.decode().strip()}")
        # Extract job ID from the output
        job_id = result.stdout.decode().strip()
        return job_id
    except Exception as e:
        print(f"Error during job submission: {e}")
        raise


def is_job_running(job_id):
    """Check if the job is still running."""
    result = subprocess.run(['qstat', job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        return True  # Job is still running
    elif "Unknown Job Id" in result.stderr:
        return False  # Job is not running (completed)
    else:
        raise RuntimeError(f"Error checking job status: {result.stderr}")


def wait_for_job_completion(job_id, poll_interval=30):
    """Wait until the job with job_id is completed."""
    while is_job_running(job_id):
        print(f"Job {job_id} is still running. Checking again in {poll_interval} seconds.")
        time.sleep(poll_interval)
    print(f"Job {job_id} has completed.")

# to load YAML as an OrderedDict
def ordered_yaml_load(file_path):
    with open(file_path, 'r') as stream:
        data = yaml.load(stream, Loader=yaml.Loader)
        return OrderedDict(data)


def prepare_directory(directory):
    """Create directory if it doesn't exist."""
    Path(directory).mkdir(parents=True, exist_ok=True)


def run_vasp_calculation(job_directory, script_path):
    """Submit VASP job and wait for completion."""
    os.chdir(job_directory)
    try:
        job_id = submit_job(script_path)
        print(f"Submitted job with ID: {job_id}")
        wait_for_job_completion(job_id)
        print("Job completed. Continuing with the next step...")
    finally:
        os.chdir(original_directory)


def setup_kpoints_hsp():
    """Setup band structure calculations."""
    structure = Structure.from_file(Path(f"{relax}/POSCAR"))
    kpath = HighSymmKpath(structure)
    kpts = Kpoints.automatic_linemode(divisions=40, ibz=kpath)
    kpts.write_file(os.path.join(f"{relax}/", "KPOINTS_hsp"))


def modify_incar_file(incar_file_path, custom_settings):
    """Modify INCAR file with custom settings."""
    incar = Incar.from_file(incar_file_path)
    #incar.update(custom_settings)
    for key, value in custom_settings.items():
        # Special handling for MAGMOM to avoid square brackets
        if key == 'MAGMOM':
            incar[key] = value.tolist()  # Convert numpy array to list
        else:
            incar[key] = value
    incar.write_file(incar_file_path)

# Copy incar files from source code
def copy_file_from_package(package_name, file_within_package, destination_path):
    # Get the file path within the package
    package_file_path = pkg_resources.resource_filename(package_name, file_within_package)
    
    # Define the destination path
    destination_file_path = Path(destination_path) / Path(file_within_package).name
    
    # Ensure the destination directory exists
    destination_file_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Copy the file to the destination
    shutil.copy(package_file_path, destination_file_path)
    print(f"Copied {package_file_path} to {destination_file_path}")

# here file_path = POTCAR file
def modify_encut(incar_file_path, file_path):
    """Modify ENCUT file with custom settings."""
    incar = Incar.from_file(incar_file_path)
    #incar.update(custom_settings)

    # Initialize a variable to store the maximum value
    max_enmax = None

    # Open the file and read its contents
    with open(file_path, 'r') as file:
        for line in file:
            # Check if the line contains 'ENMAX'
            if 'ENMAX' in line:
                # Extract the value after '=' and before ';'
                enmax_str = line.split('=')[1].split(';')[0].strip()
            
                # Convert the extracted string to float
                try:
                    enmax_value = float(enmax_str)
                except ValueError:
                    continue  # Skip this line if conversion fails
            
                # Update max_enmax if it's None or smaller than enmax_value
                if max_enmax is None or enmax_value > max_enmax:
                    max_enmax = enmax_value

    encut = int(1.5 * max_enmax)  # Adjust factor (1.3 or 1.5) as needed
    custom_settings = {"ENCUT": encut}

    for key, value in custom_settings.items():
        # Special handling for MAGMOM to avoid square brackets
        if key == 'MAGMOM':
            incar[key] = value.tolist()  # Convert numpy array to list
        else:
            incar[key] = value
    incar.write_file(incar_file_path)



# here filename will be cif_file
def extract_information(filename):
    database_code_pattern = r'^#_database_code_PCD\s+(\d+)'
    formula_pattern = r"^_chemical_formula_structural\s+'(.*?)'"

    database_code = None
    formula = None

    with open(filename, 'r') as file:
        for line in file:
            # Extract database code
            match_db_code = re.match(database_code_pattern, line)
            if match_db_code:
                database_code = match_db_code.group(1)
            
            # Extract chemical formula
            match_formula = re.match(formula_pattern, line)
            if match_formula:
                formula = match_formula.group(1)
    
    print(formula)
    # system_name = re.sub(r'[^\w\s]', '', formula)
    # system_name = system_name.replace(' ', '')
    
    # return database_code, system_name




def extract_nelect(file_path):
    nelect_value = None

    pattern = re.compile(r'NELECT\s*=\s*([\d.]+)')

    with open(file_path, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                nelect_value = int(float(match.group(1)))
                break

    return nelect_value

# file_path = './{relax}/OUTCAR'

# nelect = extract_nelect(file_path)


def modify_nbands_wsoc(incar_file_path):
    """Modify NBANDS for wsoc file with custom settings."""
    incar = Incar.from_file(incar_file_path)
    #incar.update(custom_settings)
    nelect = extract_nelect("./{relax}/OUTCAR")
    nbands = int((((int(nelect/2) + 16)/32)+1)*32)
    custom_settings = {"NBANDS": nbands}

    for key, value in custom_settings.items():
        incar[key] = value
    incar.write_file(incar_file_path)

def modify_nbands_soc(incar_file_path):
    """Modify NBANDS for soc file with custom settings."""
    incar = Incar.from_file(incar_file_path)
    #incar.update(custom_settings)
    nelect = extract_nelect("./{relax}/OUTCAR")
    nbands = int((((int(nelect) + 16)/32)+1)*32)
    custom_settings = {"NBANDS": nbands}

    for key, value in custom_settings.items():
        incar[key] = value
    incar.write_file(incar_file_path)

def modify_magmom(incar_file_path):
    """Modify MAGMOM file with custom settings."""
    incar = Incar.from_file(incar_file_path)
    #incar.update(custom_settings)
    structure = Structure.from_file(f"{soc_scf}/POSCAR")
    magmom = [0.0] * len(structure)  # Example: setting MAGMOM for all sites to 5.0 mu_B
    custom_settings = {"MAGMOM": magmom}
    incar.write_file(incar_file_path)


def plot_band_structure(vasprun_file, output_file):
    """Plot band structure."""
    vaspout = Vasprun(vasprun_file)
    bandstr = vaspout.get_band_structure(line_mode=True)
    bs_plotter = BSPlotter(bandstr)
    plot = bs_plotter.get_plot(zero_to_efermi=True, ylim=[-2, 2])
    plt.axhline(y=0, color='k', linestyle='--', linewidth=1, label="Fermi Level")
    plt.legend().set_visible(False)
    plt.savefig(output_file)
    plt.close()


def main():
    '''
    In first directory, there must exist 4 directory:
    1. cif_files
    2. calculations
    3. calculated
    4. job_files

    To change number of core there is a core_num.py
    As of now script should be inside cif_files directory. 
                or
    whenver you run the script you should be inside cif_fles directory.
    '''
    # Check if exactly one argument is provided
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <argument>")
        sys.exit(1)
    # path variables
    
    cif_file = sys.argv[1]
    
    original_directory = Path.cwd()

    #database_code, system_name = extract_information(cif_file)
    # print(f"Database code: {database_code}, System name: {system_name}")
    
    extract_information(cif_file)
    main_directory = f"{database_code}_{system_name}"
    
    relax = Path(f"../calculations/{main_directory}/relax")
    wsoc_scf = Path(f"../calculations/{main_directory}/bands/wsoc/scf")
    wsoc_dos = Path(f"../calculations/{main_directory}/bands/wsoc/dos")
    wsoc_band = Path(f"../calculations/{main_directory}/bands/wsoc/band")
    wsoc_plots = Path(f"../calculations/{main_directory}/bands/wsoc/plots")

    soc_scf = Path(f"../calculations/{main_directory}/bands/soc/scf")
    soc_dos = Path(f"../calculations/{main_directory}/bands/soc/dos")
    soc_band = Path(f"../calculations/{main_directory}/bands/soc/band")
    soc_plots = Path(f"../calculations/{main_directory}/bands/soc/plots")

    job_std = Path("../job_files/job_std.sh")
    job_ncl = Path("../job_files/job_ncl.sh")

    # Step 1: Relaxation
    struct = Structure.from_file(filename=cif_file)
    # user_incar_params = ordered_yaml_load('my_incar.yaml')
    relax_set = MPRelaxSet(structure=struct)
    relax_set.write_input(output_dir=f"{relax}")
    # Bring job_std file
    # shutil.copy(job_std, relax)
    
    #TODO
    # For incar_files
    # We will copy the INCAR file as it is from source code for the propblem 
    # of alphabetic order
    # file_to_copy = inside pymatgen source
    # target_directory = where to move
    copy_file_from_package('pymatgen', 'incar_files/INCAR_relax', f"{relax}/INCAR")

    # shutil.copy("INCAR_relax", f"{relax}/INCAR")

    # need to be changed
    modify_encut(Path(f"{relax}/INCAR"), Path(f"{relax}/POTCAR"))

    #TODO
    # Check what is ibz inside this function.
    setup_kpoints_hsp()   # Advance high symmetry path generation
    print("#####################")

    print("calculation for relax")
    run_vasp_calculation(f"{relax}", f"{job_std}")
    
    print("#####################")

    # Step 2: Band Structure SCF
    prepare_directory(f"{wsoc_scf}")
    shutil.copy(f"{relax}/CONTCAR", f"{wsoc_scf}/POSCAR")
    
    #TODO
    copy_file_from_package('pymatgen', 'incar_files/INCAR_scf', f"{wsoc_scf}/INCAR")
    
    #shutil.copy(f"{relax}/INCAR", f"{wsoc_scf}/")
    
    shutil.copy(f"{relax}/POTCAR", f"{wsoc_scf}/")
    shutil.copy(f"{relax}/KPOINTS", f"{wsoc_scf}/")
    
    # shutil.copy("job_std.sh", f"{wsoc_scf}/")

    # need to be changed
    modify_nbands_wsoc(f"{wsoc_scf}/INCAR")

    print("#####################")
    
    print("calculation for band_scf")
    run_vasp_calculation(f"{wsoc_scf}", f"{job_std}")

    print("#####################")
    

    # Step 3: Band Structure Non-SCF
    prepare_directory(f"{wsoc_band}")
    shutil.copy(f"{wsoc_scf}/POSCAR", f"{wsoc_band}/")
    shutil.copy(f"{wsoc_scf}/INCAR", f"{wsoc_band}/")
    shutil.copy(f"{wsoc_scf}/POTCAR", f"{wsoc_band}/")
    shutil.copy(f"{relax}/KPOINTS_hsp", f"{wsoc_band}/KPOINTS")
    shutil.copy(f"{wsoc_scf}/CHGCAR", f"{wsoc_band}/")
    shutil.copy(f"{wsoc_scf}/WAVECAR", f"{wsoc_band}/")
    
    # shutil.copy("job_std.sh", f"{wsoc_band}/")
    
    modify_incar_file(f"{wsoc_band}/INCAR", {'ISTART': 1,'ICHARG': 11})
    
    print("#####################")
    
    print("calculation for band_non_scf")
    run_vasp_calculation(f"{wsoc_band}", f"{job_std}")
    
    print("#####################")
    
    prepare_directory(f"{wsoc_plots}")
    plot_band_structure(f"{wsoc_band}/vasprun.xml", f"{wsoc_plots}/band.png")

    # Step 4: Band Structure with SOC_scf
    prepare_directory(soc_scf)
    shutil.copy(f"{wsoc_scf}/POSCAR", f"{soc_scf}/")
    shutil.copy(f"{wsoc_scf}/INCAR", f"{soc_scf}/")
    shutil.copy(f"{wsoc_scf}/POTCAR", f"{soc_scf}/")
    shutil.copy(f"{wsoc_scf}/KPOINTS", f"{soc_scf}/")
    
    # shutil.copy("job_ncl.sh", f"{soc_scf}/")
    
    # modify_magmom(f"{soc_scf}/INCAR")
    modify_nbands_soc(f"{soc_scf}/INCAR")
    modify_incar_file(f"{soc_scf}/INCAR", {
        'LSORBIT': '.TRUE.',
        'GGA_COMPAT': '.FALSE.',
    })

    print("#####################")
    
    print("calculation for band_soc_scf")
    run_vasp_calculation(f"{soc_scf}", f"{job_ncl}")

    print("#####################")
    

    # Step 5: Band Structure with SOC_non_scf
    prepare_directory(soc_band)
    shutil.copy(f"{soc_scf}/POSCAR", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/INCAR", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/POTCAR", f"{soc_band}/")
    shutil.copy(f"{wsoc_band}/KPOINTS", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/CHGCAR", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/WAVECAR", f"{soc_band}/")
    
    # shutil.copy("job_ncl.sh", f"{soc_band}/")
    
    modify_incar_file(f"{soc_band}/INCAR", {'ISTART': 1, 'ICHARG': 11})

    print("#####################")
    
    print("calculation for band_soc")
    run_vasp_calculation(f"{soc_band}", f"{job_ncl}")
    
    print("#####################")
    
    prepare_directory(f"{soc_plots}")
    plot_band_structure(f"{soc_band}/vasprun.xml", f"{soc_plots}/band.png")

    shutil.copy(cif_file, "calculated/")
    os.rename(f"calculated/{cif_file}", main_directory)

    print("Completed Successfully")

if __name__ == "__main__":
    main()

