# Generate DIA-ps using DIA-Umpire
# Input: Centroided DIA mzXML files in dia_data
include: "wf_diau/Snakefile"

# Generate spectral library using Philosopher
# Input: DIA-ps files in spectral_data
# Note: DDA mzXML can be copied to spectral_data and will be included
include: "wf_library/Snakefile"

# Generate quantitative matrix using OpenSWATH
# Input: Centroided DIA mzXML files in dia_data
# Input 2: Spectral library generated above
include: "wf_openswath/Snakefile"
