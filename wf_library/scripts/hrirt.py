import os
import operator
import pandas as pd

from msproteomicstoolslib.format import pepXMLReader
from sklearn import preprocessing
import statsmodels.api as sm
from scipy.interpolate import interp1d

def read_pepxml(infile):
  r = pepXMLReader.pepXMLReader(infile)

  peptides = []

  for hit in r.parse_all():
    peptides.append({'modified_peptide': hit.modified_peptide, 'rt': hit.spectrum_query.retTime, 'pp': hit.probability})

  df = pd.DataFrame(peptides)

  # Aggregate repeat identifications
  dfr = df.groupby(['modified_peptide']).apply(lambda x: x[x['pp'] == x['pp'].max()]['rt'].median()).reset_index()
  dfr.columns = ['modified_peptide','rt']

  return(dfr)

def lowess(run, reference_run):
  dfm = pd.merge(run, reference_run[['modified_peptide','irt']])

  # Fit lowess model
  lwf = sm.nonparametric.lowess(dfm['irt'], dfm['rt'], frac=.66)
  lwf_x = list(zip(*lwf))[0]
  lwf_y = list(zip(*lwf))[1]
  lwi = interp1d(lwf_x, lwf_y, bounds_error=False)

  # Apply lowess model
  run['irt'] = lwi(run['rt'])

  return run.dropna()

runs = {}
hits = {}

for pepxml in snakemake.input:
  print("Reading file %s." % pepxml)
  df = read_pepxml(pepxml)
  runs[pepxml] = df
  hits[pepxml] = df.shape[0]

# Select reference run
reference_run_id = max(hits.items(), key=operator.itemgetter(1))[0]
reference_run = runs[reference_run_id]

# Remove reference from runs
runs.pop(reference_run_id, None)

# Normalize RT of reference run
min_max_scaler = preprocessing.MinMaxScaler()
reference_run['irt'] = min_max_scaler.fit_transform(reference_run[['rt']])*100

# Normalize RT of all runs against reference
runs_irt = {}
runs_irt[reference_run_id] = reference_run
for run_id, run in runs.items():
  runs_irt[run_id] = lowess(run, reference_run)

# Write iRT calibration files
for run_id, run in runs_irt.items():
  filename = os.path.splitext(run_id)[0] + ".iRT"
  run[['modified_peptide','irt']].to_csv(filename, sep="\t", header=False, index=False)
  print("Wrote file %s containg %s peptides." % (filename, run[['modified_peptide','irt']].shape[0]))
