import os
import re
import operator
import pandas as pd

# alignment
from sklearn import preprocessing
import statsmodels.api as sm
from scipy.interpolate import interp1d

# pepXML parsing
import xml.etree.ElementTree as ET

# mzML parsing
import pymzml

def read_pepxml(infile, fdr_threshold):
  peptides = []
  namespaces = {'pepxml_ns': "http://regis-web.systemsbiology.net/pepXML"}
  ET.register_namespace('', "http://regis-web.systemsbiology.net/pepXML")
  tree = ET.parse(infile)
  root = tree.getroot()

  fdr = {}
  for error_point in root.findall('.//pepxml_ns:error_point', namespaces):
      error = float(error_point.attrib['error'])
      min_prob = float(error_point.attrib['min_prob'])
      fdr[error] = min_prob

  prob_threshold = float(fdr[min(fdr.keys(), key=lambda k: abs(k-fdr_threshold))])

  for msms_run_summary in root.findall('.//pepxml_ns:msms_run_summary', namespaces):
      base_name = os.path.basename(msms_run_summary.attrib['base_name'])
      for spectrum_query in msms_run_summary.findall('.//pepxml_ns:spectrum_query', namespaces):
        index = spectrum_query.attrib['index']
        start_scan = spectrum_query.attrib['start_scan']
        end_scan = spectrum_query.attrib['end_scan']
        assumed_charge = spectrum_query.attrib['assumed_charge']
        retention_time_sec = spectrum_query.attrib['retention_time_sec']

        for search_result in spectrum_query.findall(".//pepxml_ns:search_result", namespaces):
          for search_hit in search_result.findall(".//pepxml_ns:search_hit", namespaces):
            hit_rank = search_hit.attrib['hit_rank']

            # parse either peptide or modified peptide
            peptide = search_hit.attrib['peptide']
            protein = search_hit.attrib['protein'].split(" ")[0]
            if search_hit.find('.//pepxml_ns:modification_info', namespaces):
              modified_peptide = search_hit.find('.//pepxml_ns:modification_info', namespaces).attrib['modified_peptide']
            else:
              modified_peptide = peptide

            probability = float(search_hit.find('.//pepxml_ns:interprophet_result', namespaces).attrib['probability'])

            if probability >= prob_threshold:
              peptides.append({'base_name': base_name, 'index': int(index), 'start_scan': int(start_scan), 'end_scan': int(end_scan), 'assumed_charge': int(assumed_charge), 'retention_time_sec': float(retention_time_sec), 'hit_rank': int(hit_rank), 'peptide': peptide, 'modified_peptide': modified_peptide, 'protein': protein, 'probability': float(probability)})

  df = pd.DataFrame(peptides)
  return(df)

def lowess(run, reference_run):
  dfm = pd.merge(run, reference_run[['modified_peptide','assumed_charge','irt']])

  # Fit lowess model
  lwf = sm.nonparametric.lowess(dfm['irt'], dfm['retention_time_sec'], frac=.66)
  lwf_x = list(zip(*lwf))[0]
  lwf_y = list(zip(*lwf))[1]
  lwi = interp1d(lwf_x, lwf_y, bounds_error=False)

  # Apply lowess model
  run['irt'] = lwi(run['retention_time_sec'])

  return run

def read_mzml(mzml_path, scan_ids):
  msrun = pymzml.run.Reader(mzml_path)

  peaks_list = []
  for scan_id in scan_ids:
    spectrum = msrun[scan_id]
    peaks = pd.DataFrame(spectrum.peaks, columns=['product_mz', 'intensity'])
    peaks['precursor_mz'] = spectrum['precursors'][0]['mz']
    peaks['start_scan'] = scan_id
    peaks_list.append(peaks)

  transitions = pd.concat(peaks_list)
  return(transitions)

# Parse all pepXML files
pepid_list = []
for pepxml in snakemake.input['pepxml']:
  print("Reading file %s." % pepxml)
  pepid_file = read_pepxml(pepxml, snakemake.params['fdr_threshold'])
  pepid_list.append(pepid_file)
pepid = pd.concat(pepid_list).reset_index(drop=True)

# Patch TPP modifications
for idx, modification in pd.read_csv(snakemake.params['library_modifications']).iterrows():
  print("Replace TPP modification '%s' with UniMod modification '%s'" % (modification['TPP'], modification['UniMod']))
  pepid['modified_peptide'] = pepid['modified_peptide'].str.replace(re.escape(modification['TPP']), modification['UniMod'])

# Generate set of best replicate identifications per run
pepidr = pepid.loc[pepid.groupby(['base_name','modified_peptide','assumed_charge'])['probability'].idxmax()].sort_index()

# Select reference run
pepidr_stats = pepidr.groupby('base_name')[['modified_peptide']].count().reset_index()
print(pepidr_stats)
reference_run_base_name = pepidr_stats.loc[pepidr_stats['modified_peptide'].idxmax()]['base_name']

reference_run = pepidr[pepidr['base_name'] == reference_run_base_name]
align_runs = pepidr[pepidr['base_name'] != reference_run_base_name]

# Normalize RT of reference run
min_max_scaler = preprocessing.MinMaxScaler()
reference_run['irt'] = min_max_scaler.fit_transform(reference_run[['retention_time_sec']])*100

# Normalize RT of all runs against reference
aligned_runs = align_runs.groupby('base_name').apply(lambda x: lowess(x, reference_run)).dropna()
pepida = pd.concat([reference_run, aligned_runs]).reset_index(drop=True)

# Generate set of non-redundant global best replicate identifications
pepidb = pepida.loc[pepida.groupby(['modified_peptide','assumed_charge'])['probability'].idxmax()].sort_index()

# Prepare ID mzML pairing
peak_files = pd.DataFrame({'path': snakemake.input['mzml']})
peak_files['base_name'] = peak_files['path'].apply(lambda x: os.path.splitext(os.path.basename(x))[0])

# Parse mzML to retrieve peaks and store results in peak files
for idx, peak_file in peak_files.iterrows():
  print("Parsing file %s." % peak_file['path'])
  meta_run = pepida[pepida['base_name'] == peak_file['base_name']]
  meta_global = pepidb[pepidb['base_name'] == peak_file['base_name']]
  peaks = read_mzml(peak_file['path'], meta_run['start_scan'].tolist())
  
  # Generate run-specific PQP files
  run_pqp = pd.merge(meta_run, peaks, on='start_scan')[['precursor_mz','product_mz','intensity','irt','protein','peptide','modified_peptide','assumed_charge']]
  run_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
  run_pqp_path = os.path.splitext(peak_file['path'])[0]+"_run_peaks.tsv"
  run_pqp.to_csv(run_pqp_path, sep="\t", index=False)

  # Generate global non-redundant PQP files
  global_pqp = pd.merge(meta_global, peaks, on='start_scan')[['precursor_mz','product_mz','intensity','irt','protein','peptide','modified_peptide','assumed_charge']]
  global_pqp.columns = ['PrecursorMz','ProductMz','LibraryIntensity','NormalizedRetentionTime','ProteinId','PeptideSequence','ModifiedPeptideSequence','PrecursorCharge']
  global_pqp_path = os.path.splitext(peak_file['path'])[0]+"_global_peaks.tsv"
  global_pqp.to_csv(global_pqp_path, sep="\t", index=False)
