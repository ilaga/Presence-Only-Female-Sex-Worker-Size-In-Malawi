
library(foreign)

## Read in all 9 data sources (not spatial, those are already by cluster)
hiv.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWAR7ADT (HIV Test Results Recode)\\MWAR7AFL.DTA")

birth.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWBR7HDT (Births Recode)\\MWBR7HFL.DTA")

couple.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWCR7HDT (Couples' Recode)\\MWCR7HFL.DTA")

house.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWHR7HDT (Household Recode)\\MWHR7HFL.DTA")

ind.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWIR7HDT (Individual Recode)\\MWIR7HFL.DTA")

child.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWKR7HDT (Children's Recode)\\MWKR7HFL.DTA")

men.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWMR7HDT (Men's Recode)\\MWMR7HFL.DTA")

housemem.dat <- read.dta("D:\\Box Sync\\School\\Penn State\\Bao_Research\\DHS\\MW_2015-16_DHS\\MWPR7HDT (Household Member Recode)\\MWPR7HFL.DTA")


## Now aggregate by cluster
## The cluster variable is always "v001", except for HIV which is "hivclust"
hiv.out = c()
birth.out = c()
couple.out = c()
house.out = c()
ind.out = c()
child.out = c()
men.out = c()
housemem.out = c()

mv168 = vector(length = 850)
v155 = vector(length = 850)
v168 = vector(length = 850)
mv791 = vector(length = 850)
mv793 = vector(length = 850)
hiv.test = vector(length = 850)


for(i in 1:850){
  hiv.tmp = hiv.dat[hiv.dat$hivclust == i,]
  birth.tmp = birth.dat[birth.dat$v001 == i,]
  couple.tmp = couple.dat[couple.dat$v001 == i,]
  house.tmp = house.dat[house.dat$hv001 == i,]
  ind.tmp = ind.dat[ind.dat$v001 == i,]
  child.tmp = child.dat[child.dat$v001 == i,]
  men.tmp = men.dat[men.dat$mv001 == i,]
  housemem.tmp = housemem.dat[housemem.dat$hv001 == i,]
  
 
  ## Collect percent variables
  mv168[i] = summary(men.tmp$mv168)[1]/(summary(men.tmp$mv168)[1] + summary(men.tmp$mv168)[2])
  mv791[i] = summary(men.tmp$mv791)[1]/(summary(men.tmp$mv791)[1] + summary(men.tmp$mv791)[2])
  mv793[i] = summary(men.tmp$mv793)[1]/(summary(men.tmp$mv793)[1] + summary(men.tmp$mv793)[2])
  v155[i] = (summary(ind.tmp$v155)[1] + summary(ind.tmp$v155)[2]) /
    (summary(ind.tmp$v155)[1] + summary(ind.tmp$v155)[2] + summary(ind.tmp$v155)[3])
  v168[i] = summary(ind.tmp$v168)[1]/(summary(ind.tmp$v168)[1] + summary(ind.tmp$v168)[2])
  hiv.test[i] = summary(hiv.tmp$hiv03)[2]/(summary(hiv.tmp$hiv03)[1] + summary(hiv.tmp$hiv03)[2])
  
  print(i)
}



percent.var.out = cbind(mv168, mv791, mv793, v155, v168, hiv.test)
names(percent.var.out)[6] = "hivpos"





write.csv(hiv.out, file="hiv_agg.csv", row.names=F)
write.csv(birth.out, file="birth_agg.csv", row.names=F)
write.csv(couple.out, file="couple_agg.csv", row.names=F)
write.csv(house.out, file="household_agg.csv", row.names=F)
write.csv(ind.out, file="individual_agg.csv", row.names=F)
write.csv(child.out, file="child_agg.csv", row.names=F)
write.csv(men.out, file="men_agg.csv", row.names=F)
write.csv(housemem.out, file="housemember_agg.csv", row.names=F)
write.csv(percent.var.out, file="percent_agg.csv", row.names=F)


