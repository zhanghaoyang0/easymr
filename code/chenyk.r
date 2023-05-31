ssh compute9

library('TwoSampleMR')
library(dplyr)
df = available_outcomes()
ids = unlist(df%>%filter(grepl('bbj', id))%>%select(id))

for (id in ids[1:5]){
    url = paste0('https://gwas.mrcieu.ac.uk/files/', id, '/', id, '.vcf.gz')
    command = paste0('wget -c ', url)
    system(command, ignore.stdout=T)
}


