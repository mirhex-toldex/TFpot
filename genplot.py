from pygnuplot import gnuplot
import os
import imageio

def list_files_in_directory(directory):
    files = os.listdir(directory)
    return files

directory_path = '.'
files_in_directory = list_files_in_directory(directory_path)
sorted_files = sorted(files_in_directory)

g = gnuplot.Gnuplot()
png_file_names = []

for file in sorted_files:
    if file.startswith('fort.'):
        g.set(terminal = 'pngcairo font "arial,10" fontscale 1.0 size 600, 400', output = f'"{file}.png"')
        g.plot(f"'{file}' u 2:($3**2+$4**2) w l")
        png_file_names.append(f'{file}.png')



