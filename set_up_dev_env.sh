#!/bin/bash

if [ ! -f "set_up_dev_env.sh" ]; then
    echo "Please run this script from the root directory of the repository."
    exit 1
fi

# Path to the python virtual environment.
PYTHON_VENV="`dirname "$0"`/utilities/python-venv"

# If the virtual environment already exists, delete it.
if [ -d "$PYTHON_VENV" ]; then rm -Rf $PYTHON_VENV; fi

# Setup the virtual environment and source it.
python3 -m venv "${PYTHON_VENV}"
source "${PYTHON_VENV}"/bin/activate

# Install all the modules defined in requirements.txt.
pip install --upgrade pip
pip install wheel
pip install -r requirements.txt

# Additionally store the hash of the ingredients for the virtual environment.
./utilities/code_checks/check_venv --update

# Install the pre-commit hooks.
pre-commit install

# Copy the commit-msg hook to the .git/hooks directory.
cp utilities/code_checks/commit-msg .git/hooks/commit-msg
