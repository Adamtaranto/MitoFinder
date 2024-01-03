import shutil
import sys
import subprocess
from collections import Counter
import os.path


def is_avail(tool_names, kill=True):
    for tool_name in tool_names:
        tool_path = shutil.which(tool_name)
        if tool_path is not None:
            print(f"{tool_name} is available on the PATH: {tool_path}")
            return True
        else:
            print(f"Error: {tool_name} not found on the PATH.")
            if kill:
                sys.exit(1)
            else:
                return False


def check_files(file_names):
    for file_name in file_names:
        if file_name and os.path.exists(file_name):
            print(f"File found: '{file_name}'")
        elif file_name:
            print(f"File does not exist: '{file_name}'")
            exit(1)


def get_abs_if_found(file_name):
    if file_name and os.path.exists(file_name):
        print(f"File found: '{file_name}'")
        return os.path.abspath(file_name)
    elif file_name:
        print(f"File does not exist: '{file_name}'")
        exit(1)


def is_java_installed():
    try:
        # Run the 'java -version' command and capture the output
        java_version_output = subprocess.check_output(
            ["java", "-version"], stderr=subprocess.STDOUT, text=True
        )
        return True
    except subprocess.CalledProcessError:
        # If an error occurs, Java is likely not installed or not in the system PATH
        return False


def find_duplicates(input_list):
    # Count occurrences of each item in the list
    item_counts = Counter(input_list)

    # Filter items that occur more than once
    duplicates = [item for item, count in item_counts.items() if count > 1]

    # Remove duplicate empty strings
    duplicates = [item for item in duplicates if item]

    return duplicates


def setup_directory(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        # If it doesn't exist, create the directory
        os.makedirs(directory_path)
        print(f"Directory created: '{directory_path}'")
    else:
        # If it exists, clean existing files
        shutil.rmtree(directory_path)
        os.makedirs(directory_path)
        print(f"Directory cleaned: '{directory_path}'")
