md = read.csv("brachy-metadata.csv")
str(md)

spp = read.csv("original-metadata/species.csv", stringsAsFactors=F)
str(spp)

js = read.delim("original-metadata/Brachy_META_NO_NAs.txt")
str(js)
table(js$Plate)
table(js$Number)
table(js$Letter)

paste0(js$Letter, sprintf("%02d", js$Number))

vagrep = function(needle, haystack, ...) {
    needle = as.character(needle)
    haystack = as.character(haystack)
    res = c()
    for (n in needle) {
        r = agrep(n, haystack, ...)
        if (length(r) == 0) {
            r = NA
        }
        res = c(res, r)
    }
    res
}

m=vagrep(js$ID, md$accession, ignore.case=T)

m=vagrep(js$ID, md$accession, ignore.case=T)

