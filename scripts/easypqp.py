import sys
import os
import re
import operator
import numpy as np
import pandas as pd

# alignment
from sklearn import preprocessing
import statsmodels.api as sm
from scipy.interpolate import interp1d

# error rate estimation
from pyprophet.stats import pemp, qvalue, pi0est

# mzXML parsing
import pyopenms as po

def peptide_fdr(psms, peptide_fdr_threshold):
  pi0_lambda = np.arange(0.05, 0.5, 0.05)
  pi0_method = 'bootstrap'
  pi0_smooth_df = 3
  pi0_smooth_log_pi0 = False
  pfdr = False

  peptides = psms.groupby(['modified_peptide','decoy'])['d_score'].max().reset_index()
  targets = peptides[~peptides['decoy']]
  decoys = peptides[peptides['decoy']]

  targets['p_value'] = pemp(targets['d_score'], decoys['d_score'])
  targets['q_value'] = qvalue(targets['p_value'], pi0est(targets['p_value'], pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0)['pi0'], pfdr)
  
  return targets[targets['q_value'] < peptide_fdr_threshold]['modified_peptide']

def protein_fdr(psms, protein_fdr_threshold):
  pi0_lambda = np.arange(0.05, 0.5, 0.05)
  pi0_method = 'bootstrap'
  pi0_smooth_df = 3
  pi0_smooth_log_pi0 = False
  pfdr = False

  proteins = psms.groupby(['protein_id','decoy'])['d_score'].max().reset_index()
  targets = proteins[~proteins['decoy']]
  decoys = proteins[proteins['decoy']]

  targets['p_value'] = pemp(targets['d_score'], decoys['d_score'])
  targets['q_value'] = qvalue(targets['p_value'], pi0est(targets['p_value'], pi0_lambda, pi0_method, pi0_smooth_df, pi0_smooth_log_pi0)['pi0'], pfdr)
  
  return targets[targets['q_value'] < protein_fdr_threshold]['protein_id']

def read_psm(psm_path, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold):
  psms = pd.read_table(psm_path, index_col=False)

  # Only keep proteotypic peptides
  psms = psms[psms['proteotypic']]

  # Confident peptides and protein in global context
  peptides = peptide_fdr(psms, peptide_fdr_threshold)
  proteins = protein_fdr(psms, protein_fdr_threshold)

  # Filter peptides and proteins
  psms = psms[psms['modified_peptide'].isin(peptides)]
  psms = psms[psms['protein_id'].isin(proteins)]

  # Filter PSMs
  psms = psms[psms['q_value'] < psm_fdr_threshold]

  # Remove decoys
  psms = psms[~psms['decoy']]

  # Append columns
  psms['base_name'] = psms['run_id'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

  return psms

def lowess(run, reference_run):
  dfm = pd.merge(run, reference_run[['modified_peptide','precursor_charge','irt']])

  print("Info: Peptide overlap between run and reference: %s." % dfm.shape[0])

  # Fit lowess model
  lwf = sm.nonparametric.lowess(dfm['irt'], dfm['retention_time'], frac=.66)
  lwf_x = list(zip(*lwf))[0]
  lwf_y = list(zip(*lwf))[1]
  lwi = interp1d(lwf_x, lwf_y, bounds_error=False)

  # Apply lowess model
  run['irt'] = lwi(run['retention_time'])

  return run

def read_mzxml(mzxml_path, scan_ids):
  fh = po.MzXMLFile()
  fh.setLogType(po.LogType.CMD)
  input_map = po.MSExperiment()
  fh.load(mzxml_path, input_map)

  peaks_list = []
  for scan_id in scan_ids:

    spectrum = input_map.getSpectrum(scan_id - 1)

    product_mzs = []
    intensities = []
    for peak in spectrum:
      product_mzs.append(peak.getMZ())
      intensities.append(peak.getIntensity())

    peaks = pd.DataFrame({'product_mz': product_mzs, 'intensity': intensities})
    peaks['precursor_mz'] = spectrum.getPrecursors()[0].getMZ()
    peaks['scan_id'] = scan_id
    peaks_list.append(peaks)

  transitions = pd.concat(peaks_list)
  return(transitions)

# Parse input arguments
psms = []
mzxmls = []
psm_fdr_threshold = float(sys.argv[1])
peptide_fdr_threshold = float(sys.argv[2])
protein_fdr_threshold = float(sys.argv[3])

for arg in sys.argv[4:]:
  if 'pyprophet' in arg:
    psms.append(arg)
  if 'mzXML' in arg:
    mzxmls.append(arg)

# Parse all pepXML files
pepid_list = []
for psm in psms:
  print("Reading file %s." % psm)
  pepid_file = read_psm(psm, psm_fdr_threshold, peptide_fdr_threshold, protein_fdr_threshold)
  pepid_list.append(pepid_file)
pepid = pd.concat(pepid_list).reset_index(drop=True)

# Generate set of best replicate identifications per run
pepidr = pepid.loc[pepid.groupby(['base_name','modified_peptide','precursor_charge'])['pep'].idxmin()].sort_index()

# Select reference run
pepidr_stats = pepidr.groupby('base_name')[['modified_peptide']].count().reset_index()
print(pepidr_stats)
reference_run_base_name = pepidr_stats.loc[pepidr_stats['modified_peptide'].idxmax()]['base_name']

reference_run = pepidr[pepidr['base_name'] == reference_run_base_name]
align_runs = pepidr[pepidr['base_name'] != reference_run_base_name]

# Normalize RT of reference run
min_max_scaler = preprocessing.MinMaxScaler()
reference_run['irt'] = min_max_scaler.fit_transform(reference_run[['retention_time']])*100

# Normalize RT of all runs against reference
aligned_runs = align_runs.groupby('base_name').apply(lambda x: lowess(x, reference_run)).dropna()
pepida = pd.concat([reference_run, aligned_runs]).reset_index(drop=True)

# Generate set of non-redundant global best replicate identifications
pepidb = pepida.loc[pepida.groupby(['modified_peptide','precursor_charge'])['pep'].idxmin()].sort_index()

# Prepare ID mzML pairing
peak_files = pd.DataFrame({'path': mzxmls})
peak_files['base_name'] = peak_files['path'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

# Parse mzXML to retrieve peaks and store results in peak files
for idx, peak_file in peak_files.iterrows():
  print("Parsing file %s." % peak_file['path'])
  meta_run = pepida[pepida['base_name'] == peak_file['base_name']]
  meta_global = pepidb[pepidb['base_name'] == peak_file['base_name']]
  peaks = read_mzxml(peak_file['path'], meta_run['scan_id'].tolist())
  
  # Generate run-specific PQP files for OpenSWATH alignment
  if "_Q1" in peak_file['base_name']:
    run_pqp = pd.merge(meta_run, peaks, on='scan_id')[['precursor_mz','product_mz','intensity','irt','protein_id','peptide_sequence','modified_peptide','precursor_charge']]
    run_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
    run_pqp_path = os.path.splitext(peak_file['path'])[0]+"_run_peaks.tsv"
    run_pqp.to_csv(run_pqp_path, sep="\t", index=False)

  # Generate global non-redundant PQP files
  global_pqp = pd.merge(meta_global, peaks, on='scan_id')[['precursor_mz','product_mz','intensity','irt','protein_id','peptide_sequence','modified_peptide','precursor_charge']]
  global_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
  global_pqp_path = os.path.splitext(peak_file['path'])[0]+"_global_peaks.tsv"
  global_pqp.to_csv(global_pqp_path, sep="\t", index=False)
