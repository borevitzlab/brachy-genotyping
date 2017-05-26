library(tidyverse)

exp = read.delim("original-metadata/EXPERIMENT.tsv", stringsAsFactors=F)

brachy_exp = exp %>%
        filter(grepl("^brachy", plate.name, ignore.case=T)) %>%
        filter(grepl("^pst", BARCODE, ignore.case=T)) %>%
        filter(plate.name != "BrachyPelPlate1") %>%
        select(plate.name, BARCODE, SequencerChemistry) %>%
        mutate(plate.code = as.integer(sub('\\D*(\\d+)\\D*', '\\1', plate.name)),
               plate.code = sprintf("BR%02d", plate.code))

# Rename dupicate XT-extracted plate
brachy_exp$plate.code[brachy_exp$plate.code == "BR05" && brachy_exp$BARCODE == "PstIAdapterPlate6"] = "BR05XT"


