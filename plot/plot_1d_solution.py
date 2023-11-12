#!/usr/bin/env python3

import matplotlib.pyplot as plt
import csv
import sys

if len(sys.argv) < 2:
    print("Please enter the filename of the solution file.")
    exit(1)

datasets = []
for filename in sys.argv[1:]:
    data = None
    with open(filename) as csvfile:
        csvreader = csv.reader(csvfile)
        column_names = csvreader.__next__()
        data = dict([(col, []) for col in column_names])
        for row in csvreader:
            annotated_row = dict(zip(column_names, row))
            for k, v in data.items():
                v.append(float(annotated_row[k]))
    datasets.append(data)


def get_subplot_code(total, i):
    return total * 100 + 10 + (i + 1)

fig = plt.figure(figsize=(20,10))
color1, color2, color3 = plt.cm.viridis([0, .5, .9])
for i, data in enumerate(datasets):
    subplt = fig.add_subplot(get_subplot_code(len(datasets), i))
    ax2 = subplt.twinx()
    ax3 = subplt.twinx()

    subplt.set_title(sys.argv[1 + i])
    subplt.grid()
    subplt.set_ylabel('Density')
    ax2.set_ylabel('Velocity')
    ax3.set_ylabel('Pressure')
    subplt.set_xlabel("Element x-position [m]")

    p1 = subplt.plot(data['elemCenter'], data['rho'], label='rho', color=color1)
    p2 = ax2.plot(data['elemCenter'], data['v'], label='v', color=color2)
    p3 = ax3.plot(data['elemCenter'], data['p'], label='p', color=color3)
    subplt.legend(handles=p1 + p2 + p3, loc='best')
    
    ax3.spines['right'].set_position(('outward', 60))

plt.tight_layout()
plt.show()