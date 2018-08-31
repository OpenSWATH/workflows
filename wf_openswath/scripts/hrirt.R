library(RSQLite)
library(plyr)

input_pqp<-snakemake@input[[1]]
output_pqp<-snakemake@output[[1]]

file.copy(input_pqp, output_pqp, overwrite=TRUE)

pqp <- dbConnect(RSQLite::SQLite(), output_pqp)
anchor_candidates<-dbGetQuery(pqp, 'SELECT PRECURSOR.ID AS PRECURSOR_ID, LIBRARY_RT FROM PRECURSOR WHERE PRECURSOR.DECOY == 0;')

anchor_candidates$BIN<-cut(anchor_candidates$LIBRARY_RT, 1000, labels=FALSE)

anchors<-ddply(anchor_candidates,.(BIN),function(X){return(head(X,5))})

dbWriteTable(pqp, "temp_anchors", anchors[,c("PRECURSOR_ID","LIBRARY_RT")]) # d now in mem database

# Delete precursors
dbExecute(pqp, 'DELETE FROM PRECURSOR WHERE ID NOT IN (SELECT PRECURSOR_ID FROM temp_anchors)')

# Delete transitions
dbExecute(pqp, 'DELETE FROM TRANSITION_PRECURSOR_MAPPING WHERE PRECURSOR_ID NOT IN (SELECT PRECURSOR_ID FROM temp_anchors)')
dbExecute(pqp, 'DELETE FROM TRANSITION WHERE ID NOT IN (SELECT TRANSITION_ID FROM TRANSITION_PRECURSOR_MAPPING)')

# Delete peptides and proteins
dbExecute(pqp, 'DELETE FROM PRECURSOR_PEPTIDE_MAPPING WHERE PRECURSOR_ID NOT IN (SELECT PRECURSOR_ID FROM temp_anchors)')
dbExecute(pqp, 'DELETE FROM PEPTIDE WHERE ID NOT IN (SELECT PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING)')
dbExecute(pqp, 'DELETE FROM PEPTIDE_PROTEIN_MAPPING WHERE PEPTIDE_ID NOT IN (SELECT PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING)')
dbExecute(pqp, 'DELETE FROM PROTEIN WHERE ID NOT IN (SELECT PROTEIN_ID FROM PEPTIDE_PROTEIN_MAPPING)')

# Delete tables
dbRemoveTable(pqp,'temp_anchors')

# Clean file
dbExecute(pqp, 'VACUUM;')

dbDisconnect(pqp)
