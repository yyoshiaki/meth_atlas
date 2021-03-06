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
parser.add_argument('--reference', type=str, default='hg38', help='hg38 or mm10')
parser.add_argument('--outdir', type=str, default='./', help="output directory")
parser.add_argument('--th_cov', type=int, default=0, help='Threshold of coverage within CpG or tile.')
parser.add_argument('--bedtools', type=str, default='bedtools', help='Full path to bedtools')

args = parser.parse_args()
f = args.file
n = args.sample_name
t = args.tile
# th_cov = args.th_cov
ref = args.reference
th_cov = 5
dir_out = args.outdir
p_bedtools = args.bedtools

if not ref in ['hg38', 'mm10']:
    raise ValueError(ref + " is an invalid reference.")

def check_tile(t):
    if not os.path.exists('{d}/data/{r}.win{t}.bed.gz'.format(r=ref, d=os.path.dirname(os.path.abspath(__file__)), t=t)):
        cmd = '{b} makewindows -g {d}/data/{r}.sort.genome -w {t} > {d}/data/{r}.win{t}.bed'.format(
            b=p_bedtools, d=os.path.dirname(os.path.abspath(__file__)), t=t, r=ref)
        print(cmd)
        subprocess.run(cmd, shell=True)

        cmd = 'gzip {d}/data/{r}.win{t}.bed'.format(r=ref, d=os.path.dirname(os.path.abspath(__file__)), t=t)
        print(cmd)
        subprocess.run(cmd, shell=True)

if t != 0:
    check_tile(t)

    cmd = '{b} sort -i {f} > {f_sort}'.format(b=p_bedtools, f=f, f_sort='tmp.sort.txt')
    print(cmd)
    subprocess.run(cmd, shell=True)

    cmd = '{b} map -a {d}/data/{r}.win{x}.bed.gz -b {bis} -c 4,5,6 -o mean,sum,sum | grep -v "\.\s*\." > {o}'.format(
        b = p_bedtools, d=os.path.dirname(os.path.abspath(__file__)),
        x = t, bis='tmp.sort.txt', o='tmp.tile.txt', r=ref)
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
    # df = df[df[['meth', 'deme']].sum(axis=1) > th_cov]
    # df['methylated_frequency'] = (df['meth'] + 1) / (df['meth'] + df['deme'] + 2)
    df = df[['CpGs', 'methylated_frequency']]
    return df

df = parse_bismark(f, t)
df.columns = ['CpGs', n]
if t != 0:
    f_out = '{d}/{n}.tile{t}bp.csv'.format(d=dir_out, t=t, n=n)
else:
    f_out = '{d}/{n}.csv'.format(d=dir_out, n=n)
    df.to_csv(f_out, index=None)

df.to_csv(f_out, index=None)
print(f_out, 'generated.')

os.remove('tmp.sort.txt')
os.remove('tmp.tile.txt')
print('Temporary files deleted.')
