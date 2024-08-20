#  This source code is licensed under the MIT license. See LICENSE in the repository root directory.
#  Copyright (c) 2024. Metku team.
#  All rights reserved.

import subprocess


def upload_package():
    command = ["python", "-m", "twine", "upload", "dist/*"]
    subprocess.run(command, check=True)

if __name__ == "__main__":
    upload_package()