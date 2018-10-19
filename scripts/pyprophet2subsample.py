import pandas as pd
import sys

# Parse input arguments
psm_files = []

for arg in sys.argv[3:]:
  if 'pyprophet' in arg:
    psm_files.append(arg)

# Read all PSM files
psms_list = []
for psm_file in psm_files:
  print("Reading file %s." % psm_file)
  psms_list.append(pd.read_table(psm_file, index_col=False).sample(frac=float(sys.argv[2])/len(psm_files)))
psms = pd.concat(psms_list).reset_index(drop=True)

psms.to_csv(sys.argv[1], sep="\t", index=False)
