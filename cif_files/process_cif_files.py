import os
import subprocess
import sys

def process_cif_files(script_name):
    # Get all .cif files in the current directory
    cif_files = [f for f in os.listdir() if f.endswith('.cif')]
    
    # Sort the files to ensure consistent order
    cif_files.sort()
    
    # Check if the script exists
    if not os.path.exists(script_name):
        print(f"Error: The script '{script_name}' does not exist in the current directory.")
        return

    # Process each .cif file
    for cif_file in cif_files:
        print(f"Processing: {cif_file}")
        
        try:
            # Run the script with the .cif file as an argument
            subprocess.run([sys.executable, script_name, cif_file], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {cif_file}: {e}")
        except Exception as e:
            print(f"Unexpected error processing {cif_file}: {e}")
        
        print(f"Completed: {cif_file}\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python this_script.py <script_name>")
        sys.exit(1)

    script_name = sys.argv[1]
    process_cif_files(script_name)
