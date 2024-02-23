import os
import imageio.v2 as imageio

def list_files_in_directory(directory):
    files = os.listdir(directory)
    return files

directory_path = '.'
files_in_directory = list_files_in_directory(directory_path)
sorted_files = sorted(files_in_directory)

png_file_names = []

for file in sorted_files:
    if file.endswith('.png'):
        png_file_names.append(file)

with imageio.get_writer('animated_plot.gif', mode='I') as writer:
    for filename in png_file_names:
        image = imageio.imread(filename)
        writer.append_data(image)
        