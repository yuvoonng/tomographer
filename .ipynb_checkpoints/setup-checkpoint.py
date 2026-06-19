import subprocess
import sys
import os
from pathlib import Path

def run_setup():
    print("--- Tomographer Initialization ---")

    # 1. Install dependencies
    if os.path.exists("requirements.txt"):
        print("Installing required modules...")
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
            print("Dependencies installed successfully.")
        except subprocess.CalledProcessError:
            print("Error: Could not install dependencies. Please check your internet connection.")
            sys.exit(1)
    
    # 2. Initialize Pre-calculated Data path
    precal_data_path = Path("./Precalculated_Data")
    if precal_data_path.is_dir():
        precal_data_path = precal_data_path.resolve()
    else:
        raise ValueError("Precalculated Data does not exist.")

    line_to_add = f'fileloc = r"{precal_data_path}"\n'
    utils_file = os.path.abspath('./Runtime_Script/tomo_utils.py')
    
    if os.path.exists(utils_file):
        with open(utils_file, "r") as f:
            lines = f.readlines()
        
        if os.path.exists(utils_file):
            with open(utils_file, "r") as f:
                lines = f.readlines()
            
            # Remove any existing definition of this variable to avoid duplicates
            # even if it was previously on a different line
            lines = [line for line in lines if not line.strip().startswith("fileloc =")]
            
            # Insert the new path at index 0 (the very first line)
            lines.insert(0, line_to_add)
            
            with open(utils_file, "w") as f:
                f.writelines(lines)
    else:
        # Create the file if it doesn't exist yet
        with open(utils_file, "w") as f:
            f.write(f"{line_to_add}\n")
    
if __name__ == "__main__":
    run_setup()