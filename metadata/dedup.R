library(tidyverse)

meta = read.csv("brachy-metadata.csv") %>% 
    filter(plate.name != "BrachypodiumPlate5_PstI_XT") %>% 
    select(-original.acc.name, -Pop.name)

dups = table(meta$anon.name)
dups = names(dups)[dups > 1]

meta.good = meta %>%
    filter(!anon.name %in% dups)

meta.dup = meta %>%
    filter(anon.name %in% dups) %>% 
    unique

meta.dedup = bind_rows(meta.good, meta.dup) %>% 
    arrange(anon.name)

dups = table(meta.dedup)
dups = names(dups)[dups > 1]

write.csv(meta.dedup, "brachy-metadata.csv", row.names=F)
