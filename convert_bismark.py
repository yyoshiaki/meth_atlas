#!/usr/bin/env python

import sys
import argparse
import subprocess
import os

import pandas as pd

parser = argparse.ArgumentParser(description='convert bismark for meth_atlas')
parser.add_argument('file', type=str, help='input bismark file')
parser.add_argument('sample_name', type=str, help='sample name')
parser.add_argument('--tile', type=int, default=0, help='If tile windows, specify window size (bp).')
parser.add_argument('--bedtools', type=str, default='bedtools', help='Full path to bedtools')

args = parser.parse_args()
f = args.file
n = args.sample_name
t = args.tile
p_bedtools = args.bedtools


def check_tile(t):
    if not os.path.exists('{d}/data/hg38.win{t}.bed.gz'.format(d=os.path.dirname(os.path.abspath(__file__)), t=t)):
        cmd = '{b} makewindows -g {d}/data/hg38.sort.genome -w {t} > {d}/data/hg38.win{t}.bed'.format(
            b=p_bedtools, d=os.path.dirname(os.path.abspath(__file__)), t=t)
        print(cmd)
        subprocess.run(cmd, shell=True)

        cmd = 'gzip {d}/data/hg38.win{t}.bed'.format(d=os.path.dirname(os.path.abspath(__file__)), t=t)
        print(cmd)
        subprocess.run(cmd, shell=True)

if t != 0:
    check_tile(t)

    cmd = '{b} sort -i {f} > {f_sort}'.format(b=p_bedtools, f=f, f_sort='tmp.sort.txt')
    print(cmd)
    subprocess.run(cmd, shell=True)

    cmd = '{b} map -a {d}/data/hg38.win{x}.bed.gz -b {bis} -c 4,5,6 -o mean,sum,sum | grep -v "\.\s*\." > {o}'.format(
        b = p_bedtools, d=os.path.dirname(os.path.abspath(__file__)),
        x = t, bis='tmp.sort.txt', o='tmp.tile.txt')
    print(cmd)
    subprocess.run(cmd, shell=True)
    f = 'tmp.tile.txt'

def parse_bismark(f, t):
    df = pd.read_csv(f, sep='\t', header=None)
    df.columns = ['chromosome', 'start', 'end', 'methylated_frequency', 'meth', 'deme']
    if t == 0:
        df['CpGs'] = df['chromosome'] + ':' + (df['start']).astype(str)
    else:
        df['CpGs'] = df['chromosome'] + ':' + (df['start']).astype(str) + '-' + (df['end']).astype(str)
    df['methylated_frequency'] /= 100
    df = df[['CpGs', 'methylated_frequency']]
    return df

df = parse_bismark(f, t)
df.columns = ['CpGs', n]
if t != 0:
    f_out = '{n}.tile{t}bp.csv'.format(t=t, n=n)
else:
    f_out = '{}.csv'.format(n)
    df.to_csv(f_out, index=None)

df.to_csv(f_out, index=None)
print(f_out, 'generated.')

os.remove('tmp.sort.txt')
os.remove('tmp.tile.txt')
print('Temporary files deleted.')