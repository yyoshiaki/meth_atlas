import sys
import pandas as pd

f = sys.argv[1]
n = sys.argv[2]

def parse_bismark(f):
    df = pd.read_csv(f, sep='\t', header=None)
    df.columns = ['chromosome', 'start', 'end', 'methylated_frequency', 'meth', 'deme']
    df['CpGs'] = df['chromosome'] + ':' + (df['start']).astype(str)
    df['methylated_frequency'] /= 100
    df = df[['CpGs', 'methylated_frequency']]
    return df

df = parse_bismark(f)
df.columns = ['CpGs', n]
df.to_csv('{}.csv'.format(n), index=None)