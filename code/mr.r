cd /data/zhanghy/mr
conda activate py38_zhanghy
R


library('TwoSampleMR')
library(dplyr)
dir.create('./data/expo_dat/')
dir.create('./result/')

ao = available_outcomes()
ao = data.frame(ao)
table(ao$sex); table(ao$population)
ao1 = ao[ao$sex=='Females'&ao$population=='European',]

# save exposure_dat for each id
for (id1 in ao1$id){
    print(id1)
    file_out = paste0('./data/expo_dat/', id1, '.rdata')
    if (file.exists(file_out)){next}
    exposure_dat = extract_instruments(id1)
    save(exposure_dat, file=file_out)
}

# find valid expo
nsnp_thres = 5 # keep expo with x snp
valid_expo = c() 
for (id1 in ao1$id){
    file_in = paste0('./data/expo_dat/', id1, '.rdata')
    load(file_in)
    if (is.null(exposure_dat)){nsnp=0} else{nsnp = nrow(exposure_dat) }
    valid_expo = c(valid_expo, id1, nsnp)
}
valid_expo = data.frame(matrix(valid_expo, ncol=2, byrow=T))
names(valid_expo) = c('id', 'nsnp')
valid_expo = valid_expo[valid_expo$nsnp>=nsnp_thres,]

# mr
for (id1 in rev(valid_expo$id)){
    for (id2 in ao1$id){
        file_in = paste0('./data/expo_dat/', id1, '.rdata')
        file_out = paste0('./result/', id1, '_AND_', id2, '.rdata')
        if (id1==id2){next}
        if (file.exists(file_out)){next}
        load(file_in)
        outcome_dat = extract_outcome_data(snps=exposure_dat$SNP, outcomes = id2)
        if (is.null(outcome_dat)){res=NULL}else{
            dat = harmonise_data(exposure_dat, outcome_dat)
            res = mr(dat)
        }
        save(res, file=file_out)
    }
}




# collect mr result
nmethod_thres = 3 # keep pair with at least x method with p<0.05
collect = data.frame()
for (id1 in valid_expo$id){
    for (id2 in ao1$id){
        file_in = paste0('./result/', id1, '_AND_', id2, '.rdata')
        if (id1==id2){next}
        load(file_in)
        if (is.null(res)){next}
        res = res[res$method!='Simple mode',]
        if (sum(res$pval<0.05)>=nmethod_thres){ 
            collect = rbind(collect, res)
        }
    }
}


collect$outcome = sapply(collect$outcome, function(x){strsplit(x, ' || ', fixed=T)[[1]][1]})
collect$exposure = sapply(collect$exposure, function(x){strsplit(x, ' || ', fixed=T)[[1]][1]})
unique(paste0(collect$exposure, ' to ', collect$outcome))


collect = collect%>%filter(!(grepl('Breast cancer', exposure)&grepl('Breast cancer', outcome)))
collect = collect%>%filter(!(grepl('ovarian cancer', exposure)&grepl('ovarian cancer', outcome)))


temp = collect%>%distinct(exposure, outcome)

write.csv(collect, './temp.csv', quote=F, row.names=F)
