import shutil
import sys
import subprocess
from collections import Counter
import os.path
import logging
from typing import List


def is_avail(tool_names, kill=True):
    for tool_name in tool_names:
        tool_path = shutil.which(tool_name)
        if tool_path is not None:
            logging.info(f"{tool_name} is available on the PATH: {tool_path}")
            return True
        else:
            logging.error(f"Error: {tool_name} not found on the PATH.")
            if kill:
                sys.exit(1)
            else:
                return False


def check_files(file_names: List[str]) -> None:
    """
    Check the existence of files in a list.

    This function iterates through the list of file names and logs whether each file
    exists or not. If a file is not found, it logs an error message and exits with code 1.

    :param file_names: List of file names to check.
    :type file_names: List[str]
    :return: None
    """
    # Check if file_names is a list
    if not isinstance(file_names, list):
        raise TypeError("Input must be a list of file names.")

    for file_name in file_names:
        # Check if the file_name is not an empty string and the file exists
        if file_name and os.path.exists(file_name):
            logging.info(f"File found: '{file_name}'")
        elif file_name:
            # Log an error message and exit if the file is not found
            logging.error(f"File does not exist: '{file_name}'")
            exit(1)


def get_abs_if_found(file_name):
    if file_name and os.path.exists(file_name):
        logging.info(f"File found: '{file_name}'")
        return os.path.abspath(file_name)
    elif file_name:
        logging.error(f"File does not exist: '{file_name}'")
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
        logging.info(f"Directory created: '{directory_path}'")
    else:
        # If it exists, clean existing files
        shutil.rmtree(directory_path)
        os.makedirs(directory_path)
        logging.info(f"Directory cleaned: '{directory_path}'")


def check_if_string_in_file(f, string):
    with open(f, "r") as read_obj:
        for line in read_obj:
            if string in line:
                return True
    return False
