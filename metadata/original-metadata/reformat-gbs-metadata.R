library(tidyverse)

brachy_exp = read.csv("metadata/experiment.csv", stringsAsFactors=F) %>%
        mutate(plate.code = as.integer(sub('\\D*(\\d+)\\D*', '\\1', sample.plate.name)),
               plate.code = sprintf("BR%02d", plate.code))

samp = read.csv("metadata/sampleplatelayout.csv", stringsAsFactors=F) %>%
    filter(!is.na(plate.name)) %>% # remove the samples with missing metadata
    mutate(plate.name = sub('^Brachypodium_Plate_', 'BrachypodiumPlate', plate.name)) %>%
    filter(grepl("^brachy", plate.name, ignore.case=T)) %>%
    rename(accession=value)

# Copy of the brachy plate 5 for the XT samples
samp.xt = filter(samp, plate.name=="BrachypodiumPlate5_PstI") %>%
    mutate(plate.name="BrachypodiumPlate5_PstI_XT")
samp = rbind(samp, samp.xt)

barcodes = read.csv("metadata/barcode.csv", stringsAsFactors=F) %>%
    separate(value, c("R1Barcode", "R2Barcode"), sep='/', fill='right') %>%
    mutate(R1Barcode = sub(" ", "", R1Barcode),
           R2Barcode = sub(" ", "", R2Barcode))

brachy = samp %>%
    left_join(brachy_exp, by=c("plate.name"="sample.plate.name")) %>%
    filter(!is.na(plate.code)) %>% # remove plates not sequenced
    left_join(barcodes, by=c("barcode.plate.name"="plate.name", "well"="well")) %>%
    mutate(anon.name = paste0(plate.code, well))

writeaxe = function(DF, prefix='keyfiles/') {
    lane = unique(DF$lane)
    if (all(DF$R2Barcode == "")) {
        DF = DF[,c("R1Barcode", "anon.name")]
    } else {
        DF = DF[,c("R1Barcode", "R2Barcode", "anon.name")]
    }
    write.table(DF, paste0(prefix, lane,".axe"),
                sep='\t', quote=F, row.names=F, col.names=F)
    return(DF)
}

t = brachy %>%
    select(lane, anon.name, R1Barcode, R2Barcode) %>%
    group_by(lane) %>%
    do(writeaxe(.))

write.csv(brachy, "brachy-metadata.csv", row.names=F)
