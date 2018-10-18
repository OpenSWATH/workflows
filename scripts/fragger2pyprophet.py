import numpy as np
import pandas as pd
import sys

from Bio import SeqIO
import pyopenms as po

mdb = po.ModificationsDB()

def parse_fasta(file):
	proteins = list(SeqIO.parse(file, "fasta"))
	protein_ids = []
	protein_sequences = []
	for protein in proteins:
		protein_ids.append(protein.id)
		protein_sequences.append(str(protein.seq))
	return protein_ids, protein_sequences

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

def match_proteins(peptide, protein_ids, protein_sequences):
	decoy = True
	proteotypic = True
	indices = [i for i, s in enumerate(protein_sequences) if peptide in s]

	ids = []
	for index in indices:
		ids.append(protein_ids[index])

	return ";".join(ids)

def is_proteotypic(protein_id):
	if ";" in protein_id:
		return False
	else:
		return True

def is_decoy(protein_id):
	decoy = True
	for protein in protein_id.split(";"):
		if "DECOY" not in protein:
			decoy = False
	return decoy

fragger_names = ['scan_id','precursor_neutral_mass','retention_time','precursor_charge','rank','peptide_sequence','upstream_aa','downstream_aa','protein_id','matched_fragments','total_matched_fragments','peptide_neutral_mass','mass_difference','number_tryptic_terminii','number_missed_cleavages','modifications','hyperscore','nextscore','intercept_em','slope_em']
df = pd.read_table(sys.argv[1], header=None, names=fragger_names, index_col=False)

# Generate UniMod peptide sequence
df['modified_peptide'] = df[['peptide_sequence','modifications']].apply(match_modifications, axis=1)

# Update protein identifiers and metadata
df = df.drop(columns = 'protein_id')
protein_ids, protein_sequences = parse_fasta(sys.argv[2])
peptides = list(df['peptide_sequence'].unique())
# Find all protein identifiers
protein_ids = list(map(lambda x: match_proteins(x, protein_ids, protein_sequences), peptides))
# Assess proteotypicity
proteotypic_flags = list(map(lambda x: is_proteotypic(x), protein_ids))
# Assess decoys
decoy_flags = list(map(lambda x: is_decoy(x), protein_ids))
pepprotmap = pd.DataFrame({'peptide_sequence': peptides, 'protein_id': protein_ids, 'proteotypic': proteotypic_flags, 'decoy': decoy_flags})
df = pd.merge(df, pepprotmap, on='peptide_sequence')

# Append PyProphet columns
df['run_id'] = sys.argv[1]
df['group_id'] = df['run_id'] + "_" + df['scan_id'].astype(str)
df['expect'] = np.power(10,(df['intercept_em'] + (df['hyperscore'] * df['slope_em'])))
df['var_expectscore'] = 0.0 - np.log(df['expect'])
df['var_deltascore'] = 1.0 - (df['nextscore'] / df['hyperscore'])
df['var_lengthscore'] = np.sqrt(df['peptide_sequence'].str.len())
df = df.rename(index=str, columns={'hyperscore': 'main_var_hyperscore'})
df.to_csv(sys.argv[3], sep="\t", index=False)
