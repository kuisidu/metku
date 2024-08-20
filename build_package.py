#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2024. Metku team.
#  All rights reserved.

import subprocess
import os
import shutil

# UPDATE THIS MANUALLY
BUILD_VERSION = "0.1.16"

def build_package():
    # Run the command to build the package
    command = ['python', 'setup.py', 'sdist', "bdist_wheel"]
    subprocess.run(command, check=True)

def clean_build_artifacts():
    # Remove 'build' and 'dist' directories if they exist
    for folder in ['build', 'dist']:
        if os.path.exists(folder):
            print(f"Removing {folder} directory...")
            shutil.rmtree(folder)
            print(f"{folder} directory removed.")



if __name__ == "__main__":
    clean_build_artifacts()
    build_package()