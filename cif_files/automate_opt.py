import sys
import numpy as np
import os
import shutil
import subprocess
import time
from pathlib import Path
# import yaml
import re

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.inputs import Kpoints, Incar, Potcar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSPlotter, DosPlotter, BSPlotterProjected
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

# check what is ibz here
def setup_kpoints_hsp():
    """Setup band structure calculations."""
    structure = Structure.from_file(f"{relax}/POSCAR")
    kpath = HighSymmKpath(structure)
    kpts = Kpoints.automatic_linemode(divisions=40, ibz=kpath)
    kpts.write_file(os.path.join(f"{relax}/", "KPOINTS_hsp"))


def modify_incar_from_dict(incar_file_path, custom_settings):
    """Modify INCAR file based on custom settings dictionary."""
    # Read the existing INCAR file content
    with open(incar_file_path, 'r') as file:
        lines = file.readlines()

    # Flag to check if each key was found
    settings_found = {key: False for key in custom_settings}

    # Modify or add entries from the dictionary
    with open(incar_file_path, 'w') as file:
        for line in lines:
            # Check if the line contains any of the keys in the dictionary (even if commented)
            for key in custom_settings:
                if re.search(rf'^[!\s]*{key}\s*=', line):
                    settings_found[key] = True
                    # Comment the existing line
                    line = re.sub(rf'^[!\s]*({key}\s*=.*)', rf'!\1', line)
            file.write(line)

        # Add or update the settings based on the dictionary
        for key, value in custom_settings.items():
            if not settings_found[key]:
                file.write(f'{key} = {value}\n')
            else:
                # If the key was found, write the new value
                file.write(f'{key} = {value}\n')



# here file_path = POTCAR file
# here file_path is of POTCAR
def modify_encut(incar_file_path, file_path):
    """Modify ENCUT file with custom settings."""
    
    max_enmax = None

    with open(file_path, 'r') as file:
        for line in file:
            if 'ENMAX' in line:
                enmax_str = line.split('=')[1].split(';')[0].strip()
            
                try:
                    enmax_value = float(enmax_str)
                except ValueError:
                    continue  # Skip this line if conversion fails
            
                if max_enmax is None or enmax_value > max_enmax:
                    max_enmax = enmax_value

    if max_enmax is None:
        raise ValueError("No valid ENMAX value found in the provided file.")

    encut = int(1.5 * max_enmax)  # Adjust factor (1.3 or 1.5) as needed
    encut = ((encut + 9) // 10)*10
    # Read the existing INCAR file content
    with open(incar_file_path, 'r') as file:
        lines = file.readlines()

    # Flag to check if ENCUT was found
    encut_found = False

    # Modify or add ENCUT
    with open(incar_file_path, 'w') as file:
        for line in lines:
            # Check for ENCUT in the line (even if commented)
            if re.search(r'^[!\s]*ENCUT\s*=', line):
                encut_found = True
                # Uncomment if necessary and update the value
                line = re.sub(r'^[!\s]*ENCUT\s*=.*', f'ENCUT = {encut}', line)
            file.write(line)
        
        # If ENCUT was not found, add it to the end
        if not encut_found:
            file.write(f'\nENCUT = {encut}\n')


# here filename will be cif_file
def extract_information(filename):
    database_code_pattern = r'^#_database_code_PCD\s+(\d+)'
    formula_pattern = r"_chemical_formula_structural\s+'?(.*?)'?\s*$"

    database_code = None
    formula = None

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()  # Remove leading and trailing whitespace

            # Extract database code
            match_db_code = re.match(database_code_pattern, line)
            if match_db_code:
                database_code = match_db_code.group(1)
            
            # Extract chemical formula
            match_formula = re.search(formula_pattern, line)
            if match_formula:
                formula = match_formula.group(1)
                # If formula was found, break out of the loop
                break
    
    # If formula is still None, raise an error or handle as needed
    if formula is None:
        raise ValueError("Chemical formula not found in the file.")
    
    # Clean up the formula to create the system name
    system_name = re.sub(r'[^a-zA-Z0-9]', '', formula)
    
    # If database code is None, raise an error or handle as needed
    if database_code is None:
        raise ValueError("Database code not found in the file.")
    
    return database_code, system_name

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


def modify_nbands_wsoc(incar_file_path, outcar_file_path):
    """Modify NBANDS for wsoc file with custom settings."""
    nelect = extract_nelect(outcar_file_path)
    nbands = int(int(((int(nelect / 2) + 16) / 32) + 1) * 32)

    # Read the existing INCAR file content
    with open(incar_file_path, 'r') as file:
        lines = file.readlines()

    # Flag to check if NBANDS was found
    nbands_found = False

    # Modify or add NBANDS
    with open(incar_file_path, 'w') as file:
        for line in lines:
            # Check for NBANDS in the line (even if commented)
            if re.search(r'^[!\s]*NBANDS\s*=', line):
                nbands_found = True
                # Comment the existing NBANDS line
                line = re.sub(r'^[!\s]*(NBANDS\s*=.*)', r'!\1', line)
            file.write(line)
        
        # Add the new NBANDS value at the end
        file.write(f'NBANDS = {nbands}\n')



def modify_nbands_soc(incar_file_path, outcar_file_path):
    """Modify NBANDS for soc file with custom settings."""
    nelect = extract_nelect(outcar_file_path)
    nbands = int(int(((int(nelect) + 16) / 32) + 1) * 32)

    # Read the existing INCAR file content
    with open(incar_file_path, 'r') as file:
        lines = file.readlines()

    # Flag to check if NBANDS was found
    nbands_found = False

    # Modify or add NBANDS
    with open(incar_file_path, 'w') as file:
        for line in lines:
            # Check for NBANDS in the line (even if commented)
            if re.search(r'^[!\s]*NBANDS\s*=', line):
                nbands_found = True
                # Comment the existing NBANDS line
                line = re.sub(r'^[!\s]*(NBANDS\s*=.*)', r'!\1', line)
            file.write(line)
        
        # Add the new NBANDS value at the end
        file.write(f'NBANDS = {nbands}\n')

def dos_kpoints(input_file, output_file, reciprocal_density):
    structure = Structure.from_file(input_file)
    new_kpoints = Kpoints.automatic_density_by_vol(structure, reciprocal_density)
    new_kpoints.write_file(output_file)


# As of now not required
def modify_magmom(incar_file_path):
    """Modify MAGMOM file with custom settings."""
    incar = Incar.from_file(incar_file_path)
    #incar.update(custom_settings)
    structure = Structure.from_file(f"{soc_scf}/POSCAR")
    magmom = [0.0] * len(structure)  # Example: setting MAGMOM for all sites to 5.0 mu_B
    custom_settings = {"MAGMOM": magmom}
    incar.write_file(incar_file_path)


    relax = Path(f"../calculations/{main_directory}/relax")


def extract_and_replace_efermi_value(first_file_path, second_file_path):
    """Extracts the efermi value from the first file and replaces it in the second file."""
    
    # Extract the efermi value from the first file
    with open(first_file_path, 'r') as file:
        content1 = file.read()
    
    pattern = r'<i name="efermi">\s*(\d+\.\d+)\s*</i>'
    match = re.search(pattern, content1)
    if match:
        efermi_value = match.group(1)
    else:
        raise ValueError("efermi value not found in the first file.")
    
    # Replace the efermi value in the second file
    with open(second_file_path, 'r') as file:
        content2 = file.read()

    pattern = r'(<i name="efermi">\s*)(\d+\.\d+)(\s*</i>)'
    
    def replace_efermi(match):
        return f'{match.group(1)}{efermi_value}{match.group(3)}'

    new_content = re.sub(pattern, replace_efermi, content2)

    with open(second_file_path, 'w') as file:
        file.write(new_content)

    print(f"The efermi value {efermi_value} has been extracted from {first_file_path} and replaced in {second_file_path}.")




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

def plot_dos(vasprun_file, output_file):
    vasprun = Vasprun(vasprun_file, parse_projected_eigen=True)
    complete_dos = vasprun.complete_dos
    dos_plotter = DosPlotter()
    dos_plotter.add_dos("Total DOS", complete_dos)
    # plot = dos_plotter.get_plot()
    dos_plotter.save_plot(output_file)


def plot_projected(vasprun_file, kpoints_file, output_file):
    vasprun = Vasprun(vasprun_file, parse_projected_eigen=True)
    band_structure = vasprun.get_band_structure(kpoints_filename=kpoints_file, line_mode=True)
    bs_plotter = BSPlotterProjected(band_structure)
    # can add band_linewidth here
    bs_plotter.get_projected_plots(zero_to_efermi= True, ylim=[-2,2])
    plt.axhline(y=0, color='k', linestyle='--', linewidth=1, label="Fermi Level")
    plt.legend(loc='best')
    plt.savefig(output_file, dpi=400)

def rename_file(directory, old_filename, new_filename):
    # Construct the full file path for the old and new filenames
    old_file_path = os.path.join(directory, old_filename)
    new_file_path = os.path.join(directory, new_filename)

    # Rename the file
    os.rename(old_file_path, new_file_path)

'''
In first directory, there must exist 5 directory:
1. cif_files
2. calculations
3. calculated
4. job_files
5. incar_files

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

database_code, system_name = extract_information(cif_file)
print(f"Database code: {database_code}, System name: {system_name}")

#extract_information(cif_file)

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
def relaxation():
    struct = Structure.from_file(filename=cif_file)
    # user_incar_params = ordered_yaml_load('my_incar.yaml')
    relax_set = MPRelaxSet(structure=struct)
    relax_set.write_input(output_dir=f"{relax}")
    # Bring job_std file
    # shutil.copy(job_std, relax)

    # For incar_files
    shutil.copy("../incar_files/INCAR_relax", f"{relax}/INCAR")

    # need to be changed
    modify_encut(f"{relax}/INCAR", f"{relax}/POTCAR")
    modify_incar_from_dict(f"{relax}/INCAR", {'LWAVE': '.FALSE.'})

    #TODO
    # Check what is ibz inside this function.
    setup_kpoints_hsp()   # Advance high symmetry path generation
    print("#####################")

    print("calculation for relax")
    shutil.copy(job_std, relax)
    run_vasp_calculation(f"{relax}", "job_std.sh")

    print("#####################")

# Step 2: Band Structure SCF
def wsoc_scf_calculation():
    prepare_directory(f"{wsoc_scf}")
    shutil.copy(f"{relax}/CONTCAR", f"{wsoc_scf}/POSCAR")
    shutil.copy("../incar_files/INCAR_scf", f"{wsoc_scf}/INCAR")
    shutil.copy(f"{relax}/POTCAR", f"{wsoc_scf}/")
    shutil.copy(f"{relax}/KPOINTS", f"{wsoc_scf}/")

    # shutil.copy("job_std.sh", f"{wsoc_scf}/")

    modify_encut(f"{wsoc_scf}/INCAR", f"{wsoc_scf}/POTCAR")
    modify_nbands_wsoc(f"{wsoc_scf}/INCAR", f"{relax}/OUTCAR")

    print("#####################")

    print("calculation for band_scf")
    shutil.copy(job_std, wsoc_scf)
    run_vasp_calculation(f"{wsoc_scf}", "job_std.sh")

    print("#####################")

# Step 3: Dos calculation
def wsoc_dos_calculation():
    prepare_directory(wsoc_dos)
    shutil.copy(f"{relax}/CONTCAR", f"{wsoc_dos}/POSCAR")
    shutil.copy(f"{wsoc_scf}/INCAR", f"{wsoc_dos}/")
    shutil.copy(f"{relax}/POTCAR", f"{wsoc_dos}/")
    # shutil.copy(f"{relax}/KPOINTS", f"{wsoc_dos}/")
    dos_kpoints(f"{wsoc_dos}/POSCAR", f"{wsoc_dos}/KPOINTS", 300)

    modify_incar_from_dict(f"{wsoc_dos}/INCAR", {'ISMEAR': -5, 'LWAVE': '.FALSE.'})

    print("#####################")

    print("calculation for dos")
    shutil.copy(job_std, wsoc_dos)
    run_vasp_calculation(f"{wsoc_dos}", "job_std.sh")

    print("#####################")

    prepare_directory(wsoc_plots)
    plot_dos(f"{wsoc_dos}/vasprun.xml", f"{wsoc_plots}/dos.png")


# Step 4: Band Structure Non-SCF
def wsoc_band_calculation():
    prepare_directory(f"{wsoc_band}")
    shutil.copy(f"{wsoc_scf}/POSCAR", f"{wsoc_band}/")
    shutil.copy(f"{wsoc_scf}/INCAR", f"{wsoc_band}/")
    shutil.copy(f"{wsoc_scf}/POTCAR", f"{wsoc_band}/")
    shutil.copy(f"{relax}/KPOINTS_hsp", f"{wsoc_band}/KPOINTS")
    shutil.copy(f"{wsoc_scf}/CHGCAR", f"{wsoc_band}/")
    shutil.copy(f"{wsoc_scf}/WAVECAR", f"{wsoc_band}/")

    # shutil.copy("job_std.sh", f"{wsoc_band}/")

    modify_incar_from_dict(f"{wsoc_band}/INCAR", {'ISTART': 1,'ICHARG': 11, 'LWAVE': '.FALSE.'})

    print("#####################")

    print("calculation for band_non_scf")
    shutil.copy(job_std, wsoc_band)
    run_vasp_calculation(f"{wsoc_band}", "job_std.sh")

    print("#####################")

    extract_and_replace_efermi_value(f"{wsoc_dos}/vasprun.xml", f"{wsoc_band}/vasprun.xml")
    plot_band_structure(f"{wsoc_band}/vasprun.xml", f"{wsoc_plots}/band.png")
    # plot_projected(f"{wsoc_band}/vasprun.xml",f"{wsoc_band}/KPOINTS", f"{wsoc_plots}/band_projected.png")

# Step 5: Band Structure with SOC_scf
def soc_scf_calculation():
    prepare_directory(soc_scf)
    shutil.copy(f"{wsoc_scf}/POSCAR", f"{soc_scf}/")
    shutil.copy(f"{wsoc_scf}/INCAR", f"{soc_scf}/")
    shutil.copy(f"{wsoc_scf}/POTCAR", f"{soc_scf}/")
    shutil.copy(f"{wsoc_scf}/KPOINTS", f"{soc_scf}/")

    # shutil.copy("job_ncl.sh", f"{soc_scf}/")

    # modify_magmom(f"{soc_scf}/INCAR")
    modify_nbands_soc(f"{soc_scf}/INCAR", f"{relax}/OUTCAR")

    modify_incar_from_dict(f"{soc_scf}/INCAR", {
    'LSORBIT': '.TRUE.',
    'GGA_COMPAT': '.FALSE.',
    })

    print("#####################")

    print("calculation for band_soc_scf")
    shutil.copy(job_ncl, soc_scf)
    run_vasp_calculation(f"{soc_scf}", "job_ncl.sh")

    print("#####################")


# Step 6: Band Structure with SOC_non_scf
def soc_band_calculation():
    prepare_directory(soc_band)
    shutil.copy(f"{soc_scf}/POSCAR", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/INCAR", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/POTCAR", f"{soc_band}/")
    shutil.copy(f"{wsoc_band}/KPOINTS", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/CHGCAR", f"{soc_band}/")
    shutil.copy(f"{soc_scf}/WAVECAR", f"{soc_band}/")

    # shutil.copy("job_ncl.sh", f"{soc_band}/")

    modify_incar_from_dict(f"{soc_band}/INCAR", {'ISTART': 1, 'ICHARG': 11, 'LWAVE': '.FALSE.'})

    print("#####################")

    print("calculation for band_soc")
    shutil.copy(job_ncl, soc_band)
    run_vasp_calculation(f"{soc_band}", "job_ncl.sh")

    print("#####################")

    prepare_directory(f"{soc_plots}")
    extract_and_replace_efermi_value(f"{soc_scf}/vasprun.xml", f"{soc_band}/vasprun.xml")
    plot_band_structure(f"{soc_band}/vasprun.xml", f"{soc_plots}/band.png")


## Finally copies the completed cif file to calculated directory with same name as in calculations directory.
def final_calculation():
    shutil.copy(cif_file, "../calculated/")
    directory = "../calculated/"
    old_filename = cif_file
    new_filename = main_directory

    rename_file(directory, old_filename, new_filename)

    print("Completed Successfully")

relaxation()
wsoc_scf_calculation()
wsoc_dos_calculation()
wsoc_band_calculation()
soc_scf_calculation()
soc_band_calculation()
final_calculation()
