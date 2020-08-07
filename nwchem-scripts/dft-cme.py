import csv
import os
import sys
import argparse
from multiprocessing import cpu_count
from os.path import *
import pandas as pd
from pkg_resources import resource_filename
from snakemake import snakemake

#
#   Create initial file structure
#

def csv2inchi(inchilist):

    with open(inchilist,'r') as f:
        if 'test.csv' in inchilist:
            reader = csv.reader(f,delimiter=(','))

        lineCount = 0
        for row in reader:
            if lineCount == 0:
                InChIKey = row.index("InChI-Key")
                InChI = row.index("InChI")
                lineCount+=1
            else:
                print(row)
                inchikey_str=row[InChIKey].split('=')[1]
                print('inchikey string ', inchikey_str)

                moldir = inchikey_str
                if not os.path.exists(moldir):
                    os.mkdir(moldir)


                initial_structure_dir = moldir + '/initial_structure'
                if not os.path.exists(initial_structure_dir):
                        os.mkdir(initial_structure_dir)

                md_structure_dir = moldir + '/md'
                if not os.path.exists(md_structure_dir):
                        os.mkdir(md_structure_dir)

                dft_structure_dir = moldir + '/dft'
                if not os.path.exists(dft_structure_dir):
                        os.mkdir(dft_structure_dir)

                inchifile_str = initial_structure_dir + '/' + inchikey_str + '.inchi'
                with open(inchifile_str,'w+') as f:
                    f.write(row[InChI])

csv2inchi('test.csv')
os.system('snakemake -p --cores 2 --snakefile snakemake/final_pipeline.snakemake -w 300')
