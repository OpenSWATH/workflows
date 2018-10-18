import numpy as np
import pandas as pd
import sys

import pyopenms as po

mdb = po.ModificationsDB()

def match_modifications(peptide):
	modified_peptide = peptide['peptide_sequence']
	modifications = {}
	if "M|" in peptide['modifications']:
		for modification in peptide['modifications'].split('|')[1:]:
			site, mass = modification.split('$')
			modifications[int(site)] = float(mass)

		for site in sorted(modifications, reverse=True):
			rm = mdb.getBestModificationByDiffMonoMass(modifications[site], 0.02, peptide['peptide_sequence'][site], po.ResidueModification.TermSpecificity.ANYWHERE)
			modified_peptide = modified_peptide[:site+1] + "(UniMod:" + str(rm.getUniModRecordId()) + ")" + modified_peptide[site+1:]

	return modified_peptide

fragger_names = ['scan_id','precursor_neutral_mass','retention_time','precursor_charge','rank','peptide_sequence','upstream_aa','downstream_aa','protein_id','matched_fragments','total_matched_fragments','peptide_neutral_mass','mass_difference','number_tryptic_terminii','number_missed_cleavages','modifications','hyperscore','nextscore','intercept_em','slope_em']
df = pd.read_table(sys.argv[1], header=None, names=fragger_names, index_col=False)

# Generate UniMod peptide sequence
print("Info: Matching modifications to UniMod.")
df['modified_peptide'] = df[['peptide_sequence','modifications']].apply(match_modifications, axis=1)

# Update protein identifiers and metadata
print("Info: Matching peptides to proteins.")
df = df.drop(columns = 'protein_id')
peptide_index = pd.read_pickle(sys.argv[2])
peptides_1 = df.shape[0]
df = pd.merge(df, peptide_index, on='peptide_sequence')
peptides_2 = df.shape[0]

if peptides_1 != peptides_2:
	sys.exit("Error: Peptides from PSMs (%s) don't match peptides after matching with FASTA (%s). Check digestion parameters." % (peptides_1, peptides_2))

# Append PyProphet columns
df['run_id'] = sys.argv[1]
df['group_id'] = df['run_id'] + "_" + df['scan_id'].astype(str)
df['expect'] = np.power(10,(df['intercept_em'] + (df['hyperscore'] * df['slope_em'])))
df['var_expectscore'] = 0.0 - np.log(df['expect'])
df['var_deltascore'] = 1.0 - (df['nextscore'] / df['hyperscore'])
df['var_lengthscore'] = np.sqrt(df['peptide_sequence'].str.len())
df = df.rename(index=str, columns={'hyperscore': 'main_var_hyperscore'})
df.to_csv(sys.argv[3], sep="\t", index=False)
