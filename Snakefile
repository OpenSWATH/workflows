# Generate DIA-ps using DIA-Umpire
# Input: Centroided DIA mzXML files in dia_data
subworkflow diau:
    workdir: "wf_diau"

# Generate spectral library using Philosopher
# Input: DIA-ps files in spectral_data
# Note: DDA mzXML can be copied to spectral_data and will be included
subworkflow library:
    workdir: "wf_library"

# Generate quantitative matrix using OpenSWATH
# Input: Centroided DIA mzXML files in dia_data
# Input 2: Spectral library generated above
subworkflow openswath:
    workdir: "wf_openswath"

rule all:
    input: 
        diau = diau("wf_diau.done"),
        library= library("wf_library.done"),
        openswath = openswath("wf_openswath.done")
