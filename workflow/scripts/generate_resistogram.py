# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

# import required packages
import pandas as pd
import numpy as np
import re
import argparse
import textwrap

class Formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.MetavarTypeHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

parser = argparse.ArgumentParser(formatter_class=Formatter, 
                                description=textwrap.dedent('''
         Please do not mess up this text!
         --------------------------------
        Script for estimating efficency of antibodies against SARS-CoV-2 infection in the presence of mutations. 
        If the patient number is not in the list all known antibodies are efficient. 
        Predictions are based on the work of Lui et.al (https://doi.org/10.1038/s41586-021-04388-0), who investigated the impact of individual mutations within B.1.1.529 against monoclonal antibodies. 
        Laboratory data regarding mutation interactions and interdependencies is unavailable. Therefore, the presented algorithm is not able to interfere possible mutation interactions and interdependencies from the data.
        Since mutation interactions are expected the suggested antibody might be wrong.'''), epilog='Please cite me if you are using this program. Have fun! :)')

parser.add_argument('-p', "--print_string", help="Prints the supplied argument.", default = 'A random string.', nargs='?', const='A random string.', type=str) #works with value after flag given or not 
parser.add_argument('-k', "--keks", help="Prints the supplied argument.", type=int) #expects value after flag
parser.add_argument('-z', "--zoom", help="Prints the supplied argument.", action="store_true") #flag sets boolean to True
parser.add_argument('-w', '--week', help="Selects timespan", default = 'week87', type=str)
parser.add_argument('-f', '--frequency', help="The frequency, or allel frequency, is the relative frequency of a variant of a nucleic acid at a particular position in a genome.", default = 0.0, nargs='?', const=0.0, type=float)
parser.add_argument('-d', '--readdepth', help='Number of reads at a particular position in a genome', default = 10, nargs='?', const = 10, type=int)
args = parser.parse_args()

#retrive args
print(args)
week = str(args.week)
freq = float(args.frequency)
rd = int(args.readdepth)

# list of mutations effectiv in escaping the antibodies.
escaping_mutations = pd.read_json('/homes/kblock/scripts/resistogram/mabs.json').set_index('mAbs')
all_escaping_mutations = [item for sublist in escaping_mutations['Mutation'].tolist() for item in sublist]

# list of factors
factors = pd.read_csv("/groups/ds/kblock/resistogram/factortab.csv")
factors = factors[['Mutation', 'S309', 'COV2-2130', 'COV2-2196']].set_index('Mutation')

# list of mutations in samples
df = pd.read_csv("/groups/ds/kblock/virusrecombination/results/allmutationsfound.csv")
# filter for tests
# TODO keep all filenames ?list?
df = df[(df['Gen'] == 'S') & (df['Week'] == week) & (df['Frequency'] > freq) & (df['ReadDepth'] > rd) & ((df['ReadDepth'] * df['Frequency']) > 1)]
# filter of all samples containing concerning mutations
df = df[df['Signature'].str.contains('|'.join(all_escaping_mutations),case=False, na = False)].reset_index(drop=True)
df = df.drop(columns= ['Position', 'ReadDepth', 'Probability'])
df = df[['Filename', 'Signature', 'Gen', 'Week', 'Frequency', 'Lineage']].reset_index(drop=True)
df["ShortSig"] = np.nan
df["Ineffektiv"] = np.nan

# lists containing mutations 
medis = escaping_mutations.index.values
S309 = np.array([escaping_mutations["Mutation"]["S309"]]) #["S:I332","S:T333","S:N334","S:L335","S:C336","S:P337","S:G339","S:E340","S:V341","S:N343","S:A344","S:T345","S:R346","S:N354","S:K356","S:R357","S:I358","S:S359","S:N360","S:C361","S:N440","S:L441","S:R509"]
AZD1061 = np.array([escaping_mutations["Mutation"]["AZD1061"]]) #["S:T345","S:R346","S:N439","S:N440","S:L441","S:S443","S:K444","S:V445","S:G446","S:G447","S:Y449","S:N450","S:E484","S:F490","S:L492","S:Q493","S:S494","S:P499"]
AZD8895 = np.array([escaping_mutations["Mutation"]["AZD8895"]]) #["S:L455","S:F456","S:A475","S:G476","S:S477","S:T478","S:P479","S:E484","S:G485","S:F486","S:N487","S:C488","S:Y489","S:Q493"]

#create data structures for loop
df_sid = pd.DataFrame(columns=['Sample_id', 'Antibody', 'hig_imp_fac', 'Escaping_mutations'])
df_g = pd.DataFrame(columns=['Sample_id', 'Antibody', 'hig_imp_fac', 'Escaping_mutations'])
shortsigs =[]

#creation of output file
for idx, (sid, group) in enumerate(df.groupby('Filename')):
    #check if there are any antibodies against which there are no mutations? <- Check with Folker
    group = group.reset_index(drop=True)
    for ix, sig in enumerate(group['Signature']):
        short_sig = re.split('(\d+)', sig)[0] + re.split('(\d+)', sig)[1]
        shortsigs.append(short_sig)
    i_fac = factors[factors.index.str.contains('|'.join(shortsigs),case = False, na = False)]
    n_fac = i_fac.applymap(lambda x: (x - (-1000))/(8.6 - (-1000)))
    ab = n_fac.sum().sort_values(ascending = False)
    for i, a in enumerate(ab.index):
        ls = sorted(list(zip(i_fac[ab.index[i]].index, i_fac[ab.index[i]])), key=lambda a: a[1], reverse = True)
        df_g.at[i, "Escaping_mutations"] = ls
        df_g.at[i, "Antibody"] = a
        df_g.at[i, "Sample_id"] = sid
    df_sid = pd.concat([df_sid, df_g])

print(df_sid)
df_sid.to_csv('resistogram.csv', index=False)