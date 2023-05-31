
# About easymr
Find potential causal relationships between a numerious  exposures and outcomes is a command need.  
`easymr` is a simple pipeline to do so.

# Requirements 
- `Linux` 
- `R` (>=4.0) with `TwoSampleMR`, `dplyr` 

# Pipeline
Clone this repository via the commands:
```  
git clone https://github.com/zhanghaoyang0/easymr.git
cd easymr
R
```

Load packages and functions:
```
library('TwoSampleMR')
library(dplyr)
is_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}
```

Below is a pipeline to perform MR analysis on neonatal diseases and diabetes.

Select exposures and outcomes:
```
# load all traits
ao = data.frame(available_outcomes())
sprintf('number of all traits: %.0f', nrow(ao))
print('frq of sex:')
table(ao$sex)
# Females             Males Males and Females                NA 
# 242              3185             36335              2584 
print('frq of population:')
table(ao$population)
# European       South Asian       ...
# 38065       1503       ...

# select population
selected_sex= 'Males and Females'
selected_population = 'European'

ao1 = ao%>%filter(sex==selected_sex&population==selected_population)
sprintf('number of traits in %s sex and in %s population: %.0f', selected_sex, selected_population, nrow(ao1))
# number of traits in Males and Females sex and in European population: 32336

# select expo traits and outcome traits by patterns (keyword)
selected_expo = 'Neona|child|Child' # expo name should contain either Neona, child, Child
selected_outcome = 'Type 2 diabetes' # outcome name should contain either diabetes, Diabetes

expo_id = ao1%>%filter(grepl(selected_expo, trait))%>%pull(id)
outcome_id = ao1%>%filter(grepl(selected_outcome, trait))%>%pull(id)
outcome_id = setdiff(outcome_id, expo_id)
sprintf('number of exposure traits: %.0f', length(expo_id))
sprintf('number of outcome traits: %.0f', length(outcome_id))
# number of exposure traits: 42
# number of outcome traits: 15
```

Get exposure data:
```
for (id1 in expo_id){
    print(paste0('extracting expo_dat for ', id1))
    file_out = paste0('./data/expo_dat/', id1, '.rdata')
    if (file.exists(file_out)){next}
    exposure_dat = extract_instruments(id1)
    save(exposure_dat, file=file_out)
}

nsnp_thres = 10 # keep expo with x snp
out = c() 
for (id1 in expo_id){
    file_in = paste0('./data/expo_dat/', id1, '.rdata')
    load(file_in)
    if (is.null(exposure_dat)){nsnp=0} else{nsnp = nrow(exposure_dat)}
    out = c(out, id1, nsnp)
}
res = data.frame(matrix(out, ncol=2, byrow=T))
res = res%>%rename('id'=X1, 'nsnp'=X2)%>%mutate(nsnp=as.numeric(nsnp))
valid_expo_id = res%>%filter(nsnp>=nsnp_thres)%>%pull(id)
sprintf('number of valid exposure traits: %.0f', length(valid_expo_id))
# number of valid exposure traits: 3
```


Perform MR:
```
for (id1 in valid_expo_id){
    for (id2 in outcome_id){
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
methods = c('IVW', 'MrEgger', 'WeightedMedian', 'WeightedMode')
out = c()
for (id1 in valid_expo_id){
    for (id2 in outcome_id){
        file_in = paste0('./result/', id1, '_AND_', id2, '.rdata')
        if (id1==id2){next}
        load(file_in)
        if (is.null(res)){next}
        res = res[res$method!='Simple mode',]
        res[res=='Inverse variance weighted'] = 'IVW'
        res[res=='MR Egger'] = 'MrEgger'
        res[res=='Weighted median'] = 'WeightedMedian'
        res[res=='Weighted mode'] = 'WeightedMode'
        out = c(out, res$exposure[1], res$outcome[2])
        for (method_ in methods){
            stat = unlist(res%>%filter(method==method_)%>%select(nsnp, b, se, pval))
            out = c(out, stat)
        }
    }
}

res = data.frame(matrix(out, ncol=18, byrow=T))%>%mutate_if(is_numeric,as.numeric)
names(res) = c('exposure', 'outcome', paste0(rep(methods, each=4), rep(c('_nsnp', '_beta', '_se', '_p'), length(methods))))

# you may pick some to see:
res%>%filter(IVW_p<0.05)

#            exposure                                       outcome   IVW_nsnp  IVW_beta    IVW_se      IVW_p MrEgger_nsnp MrEgger_beta MrEgger_se ...
# Childhood sunburn occasions || id:ukb-b-13246 Type 2 diabetes || id:ieu-a-26       47 -0.447458 0.1996086 0.02498229           47    0.2774541  0.3745701 ...
```


# Feedback and comments
Add an issue or send email to zhanghaoyang0@hotmail.com