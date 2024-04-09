# STEP 4 - Create the data to be imported to theStringApp in Cytoscape
# -----------------------------------------------------------------------------
#
# This script will create and write the data needed to be imported into Cytoscape, and used in the StringApp in order to create a PPI Network. And then to perform clustering on the network.
#
# The file is stored in `data/derived-data/`. 

# Create the data
DR_output4string <- merge(BMDres_definedCI, dr_t_regs, 
                          by.x = "id", by.y = "ensembl_transcript_id_version")

# Save the data
write.table(DR_output4string, file = "outputs/DR_output4string.txt", row.names = FALSE, sep = "\t")

# Once the clustered network is created, the resulting *.csv* files need to be stored in `data/derived-data/`.
