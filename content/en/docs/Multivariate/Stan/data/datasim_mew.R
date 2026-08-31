setwd("C:/Dropbox/WAMBAM/hackathon 2026/mews")

#sim pedigree
library(pedAgree)
ped = pedAgree::simulate_pedigree(years = 10)
ped = ped$pedigree
write.csv(ped, "mewped.csv", row.names = F)

#sim phenotypes
library(mvtnorm)
ids = ped[ped$cohort > 7, "animal"]
lid = length(ids)
phens = data.frame(animal = ids, sex = sample(c(0,1), 
                   lid, replace = T))
means = c(5, 8, 30) #body mass (kg), tail length (cm), shyness (sec)
vars = (means * 0.2)^2
gvars = vars * c(0.5, 0.7, 0.3)
evars = vars - gvars
gsds = diag(sqrt(gvars))
esds = diag(sqrt(evars))
gcor = matrix(c(1, 0.6, -0.5,
               0.6, 1, -0.4,  #size-based plasticity 
               -0.5, -0.4, 1), 3, 3) #+size leads to -shyness
ecor = matrix(c(1, 0.4, 0.2, #greater integration in
                0.4, 1, 0.1, #+ resource environments
                0.2, 0.1, 1), 3, 3)
gcov = gsds %*% gcor %*% gsds
ecov = esds %*% ecor %*% esds

#genetic and environmental values
library(nadiv)
A = as.matrix(makeA(ped[,c("animal","dam","sire")]))
A = A[phens$animal, phens$animal] #subset to match phenotypic data
astd = rmvnorm(lid, mean=rep(0,3), sigma=diag(rep(1,3)))
a = t(chol(A)) %*% astd %*% t( gsds %*% t(chol(gcor)) )
cor(a)

e = rmvnorm(lid, mean = rep(0,3), sigma = ecov)
cor(a)
cor(e)
p = a + e

Beta = matrix(c(0.5, 0, -5), 1, 3)
p = a + e
p = sweep(p, 2, means, "+")
#to add fixed effects of interest, e.g. sex
#p = p + as.matrix(phens[,c("sex")]) %*% Beta
apply(p, 2, range) #sensible estimates

phens$mass_kg = p[,1]
phens$tail_cm = p[,2]
phens$shy_sec = p[,3]
phens$sex = ifelse(phens$sex == 0, "F", "M")

#save phenotypes
write.csv(phens, "mews.csv", row.names = F)

#save true pop and sample values
tdf = list(Gv = gvars, Gr = gcor, G = gcov, 
           Ev = evars, Er = ecor, E = ecov,
           h2 = gvars / (gvars + evars),
           B = Beta, a = a, e = e)
saveRDS(tdf, "mews_simval.rds")
