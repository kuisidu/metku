#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2024. Metku team.
#  All rights reserved.

import subprocess
import os
import shutil
import sys

# UPDATE THIS MANUALLY
BUILD_VERSION = "0.1.21"


def run_command(command):
    """Run a shell command and return the output."""
    try:
        result = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}\n{e.stderr}")
        sys.exit(1)


def check_merge_conflicts():
    """Check if there are merge conflicts."""
    conflicts = run_command(['git', 'ls-files', '-u'])
    if conflicts:
        print("Merge conflicts detected. Please resolve them manually.")
        sys.exit(1)


def check_uncommitted_changes():
    """Check if there are any uncommitted changes."""
    status = run_command(['git', 'status', '--porcelain'])
    return bool(status)


def stash_changes():
    """Stash uncommitted changes."""
    print("Stashing uncommitted changes...")
    run_command(['git', 'stash'])


def apply_stash():
    """Apply stashed changes if available."""
    stash_list = run_command(['git', 'stash', 'list'])
    if stash_list:
        print("Applying stashed changes...")
        run_command(['git', 'stash', 'pop'])
    else:
        print("No stashed changes to apply.")


def build_package():
    # Run the command to build the package
    # command = ['python', 'setup.py', 'sdist', "bdist_wheel"]
    command = ['python', '-m', 'build', '--sdist', '--wheel']
    subprocess.run(command, check=True)

def clean_build_artifacts():
    # Remove 'build' and 'dist' directories if they exist
    for folder in ['build', 'dist']:
        if os.path.exists(folder):
            print(f"\rRemoving {folder} directory...")
            shutil.rmtree(folder)
            print(f"\r{folder} directory removed.")


def main():
    # Step 1: Pull the newest changes
    print("Pulling the latest changes...")
    run_command(['git', 'pull', 'origin', 'master'])

    # Step 2: Check for merge conflicts
    check_merge_conflicts()

    # Step 3: Stash uncommitted changes if they exist
    # TODO: Doesnt work properly
    # if check_uncommitted_changes():
    #     stash_changes()
    # else:
    #     print("No uncommitted changes to stash.")

    # Step 4: clean build artifacts
    clean_build_artifacts()

    # Step 5: Build the repository
    build_package()

    # Step 6: Retrieve stashed changes
    # TODO: Doesnt work properly
    # apply_stash()

    print("Pipeline completed successfully!")

if __name__ == "__main__":
    main()