library(LFQbench)

sampleComposition = data.frame( 
  species = c("HUMAN","YEAST", "ECOLI"), 
  A       = c(  67,     30,       3   ), 
  B       = c(  67,      3,      30   )
)

dataSets = data.frame(
  "HYE_LUMOS" = c(
    "samon_D1701_167", "samon_D1701_168", "samon_D1701_169", # A
    "samon_D1701_172", "samon_D1701_173", "samon_D1701_175"  # B
  ),
  row.names = c( "A1", "A2", "A3", "B1", "B2", "B3" )
)

speciesTags = list(
  HUMAN = "_HUMAN", 
  YEAST = "_YEAST", 
  ECOLI = "_ECOLI"
)

LFQbench.initConfiguration(
  SampleComposition = sampleComposition
)

FSWE.initConfiguration( 
  injectionNames = dataSets,
  speciesTags = speciesTags
)

srcDir = dirname(snakemake@input[[1]])

LFQbench.setDataRootFolder( 
  rootFolder = srcDir, 
  createSubfolders = T
)

inputFiles = list.files(
  path = LFQbench.Config$DataRootFolder, 
  pattern = "\\..+"
)

nix = sapply(
  inputFiles, 
  FSWE.generateReports,
  softwareSource = "guess",
  keep_original_names = T,
  singleHits = F,
  plotHistogram = T, 
  plotHistNAs = T, 
  reportSequences = F
)

hye.res = LFQbench.batchProcessRootFolder()
