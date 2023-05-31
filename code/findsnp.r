library('TwoSampleMR')
library(dplyr)

ao = data.frame(available_outcomes())
ao = ao%>%select(id, trait, ncase, ncontrol, year, author, population, sex, pmid)
ao$snp = NA

pthres = 5e-8
for (i in 1:nrow(ao)[1:5]){
    print(i)
    id = ao[i, 'id']
    snp = extract_instruments(id, p1 = pthres, clump = F)
    snp = paste0(snp$SNP, collapse=';')
    ao[i, 'snp'] = snp
}

ao