import os

def list_files_in_directory(directory):
    files = os.listdir(directory)
    return files

directory_path = '.'
files_in_directory = list_files_in_directory(directory_path)

for file in files_in_directory:
    if file.startswith("fort."):
        os.remove(file)