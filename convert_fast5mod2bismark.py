#!/usr/bin/env python

import sys
import argparse
import subprocess
import os

import pandas as pd

parser = argparse.ArgumentParser(description='convert fast5mod output to bismark.cov')
parser.add_argument('file', type=str, help='input fast5mod output such as meth.tsv')
parser.add_argument('output', type=str, help='output file. ex) Fr2.bismark.cov.gz')

args = parser.parse_args()
f_i = args.file
f_o = args.output

df_nanopore = pd.read_csv(f_i, sep='\t', header=None)
df_nanopore.columns = ['chr', 'pos', 'CG', 'fwd_meth', 'rev_meth', 'fwd_unmeth', 'rev_unmeth']
df_nanopore['pos'] += 1
df_nanopore['meth'] = df_nanopore['fwd_meth'] + df_nanopore['rev_meth']
df_nanopore['unmeth'] = df_nanopore['fwd_unmeth'] + df_nanopore['rev_unmeth']
df_nanopore = df_nanopore[df_nanopore['meth'] + df_nanopore['unmeth'] > 0]
df_nanopore['rate'] = df_nanopore['meth'] / (df_nanopore['meth'] + df_nanopore['unmeth'])
df_nanopore = df_nanopore[['chr', 'pos', 'pos', 'rate', 'meth', 'unmeth']]
df_nanopore.columns = ['chr', 'pos', 'pos2', 'rate', 'meth', 'unmeth']
df_nanopore.to_csv(f_o, sep='\t', header=None, index=None)