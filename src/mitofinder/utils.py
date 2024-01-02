import shutil
import sys


def is_avail(tool_names):
    for tool_name in tool_names:
        tool_path = shutil.which(tool_name)
        if tool_path is not None:
            print(f"{tool_name} is available on the PATH: {tool_path}")
        else:
            print(f"Error: {tool_name} not found on the PATH.")
            sys.exit(1)
