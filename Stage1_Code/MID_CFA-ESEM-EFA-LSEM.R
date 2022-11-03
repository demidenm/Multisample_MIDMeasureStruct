# Automatically check and install necessary packages that were used in R v 4.0
if (!require("pacman")) install.packages("pacman")
pacman::p_load(simsem, plyr, tidyverse,corrplot, reshape2, data.table,ggplot2,esemComp,
               Hmisc,nFactors,car,psych, paran, #for EFA
               semTools, semPlot, lavaan, # FOR CFA/ESEM
               parameters, superheat,weights, devtools, heatmaply, # For plotting/tabling
               sirt, msm)

# Simulating data
## Specify Population Model
population_model<-'
# By run loadings for bilateral regions
AWin_v_Neut_L_NAcc =~     .7*AWin_v_Neut_L_NAcc_run1    + .7*AWin_v_Neut_L_NAcc_run2
AWin_v_Neut_L_Insula =~   .7*AWin_v_Neut_L_Insula_run1  + .7*AWin_v_Neut_L_Insula_run2
BWin_v_Neut_L_NAcc =~     .7*BWin_v_Neut_L_NAcc_run1    + .7*BWin_v_Neut_L_NAcc_run2
BWin_v_Neut_L_Insula =~   .7*BWin_v_Neut_L_Insula_run1  + .7*BWin_v_Neut_L_Insula_run2
BWin_v_BLose_L_NAcc =~    .7*BWin_v_BLose_L_NAcc_run1   + .7*BWin_v_BLose_L_NAcc_run2
BWin_v_BLose_L_Insula =~  .7*BWin_v_BLose_L_Insula_run1 + .7*BWin_v_BLose_L_Insula_run2
ALose_v_Neut_L_NAcc =~    .7*ALose_v_Neut_L_NAcc_run1   + .7*ALose_v_Neut_L_NAcc_run2
ALose_v_Neut_L_Insula =~  .7*ALose_v_Neut_L_Insula_run1 + .7*ALose_v_Neut_L_Insula_run2
BLose_v_Neut_L_NAcc =~    .7*BLose_v_Neut_L_NAcc_run1   + .7*BLose_v_Neut_L_NAcc_run2
BLose_v_Neut_L_Insula =~  .7*BLose_v_Neut_L_Insula_run1 + .7*BLose_v_Neut_L_Insula_run2
BLose_v_BWin_L_NAcc =~    .7*BLose_v_BWin_L_NAcc_run1   + .7*BLose_v_BWin_L_NAcc_run2
BLose_v_BWin_L_Insula =~  .7*BLose_v_BWin_L_Insula_run1 + .7*BLose_v_BWin_L_Insula_run2

AWin_v_Neut_R_NAcc =~     .7*AWin_v_Neut_R_NAcc_run1    + .7*AWin_v_Neut_R_NAcc_run2
AWin_v_Neut_R_Insula =~   .7*AWin_v_Neut_R_Insula_run1  + .7*AWin_v_Neut_R_Insula_run2
BWin_v_Neut_R_NAcc =~     .7*BWin_v_Neut_R_NAcc_run1    + .7*BWin_v_Neut_R_NAcc_run2
BWin_v_Neut_R_Insula =~   .7*BWin_v_Neut_R_Insula_run1  + .7*BWin_v_Neut_R_Insula_run2
BWin_v_BLose_R_NAcc =~    .7*BWin_v_BLose_R_NAcc_run1   + .7*BWin_v_BLose_R_NAcc_run2
BWin_v_BLose_R_Insula =~  .7*BWin_v_BLose_R_Insula_run1 + .7*BWin_v_BLose_R_Insula_run2
ALose_v_Neut_R_NAcc =~    .7*ALose_v_Neut_R_NAcc_run1   + .7*ALose_v_Neut_R_NAcc_run2
ALose_v_Neut_R_Insula =~  .7*ALose_v_Neut_R_Insula_run1 + .7*ALose_v_Neut_R_Insula_run2
BLose_v_Neut_R_NAcc =~    .7*BLose_v_Neut_R_NAcc_run1   + .7*BLose_v_Neut_R_NAcc_run2
BLose_v_Neut_R_Insula =~  .7*BLose_v_Neut_R_Insula_run1 + .7*BLose_v_Neut_R_Insula_run2
BLose_v_BWin_R_NAcc =~    .7*BLose_v_BWin_R_NAcc_run1   + .7*BLose_v_BWin_R_NAcc_run2
BLose_v_BWin_R_Insula =~  .7*BLose_v_BWin_R_Insula_run1 + .7*BLose_v_BWin_R_Insula_run2

#Factor item loadings 
Approach =~  .8*AWin_v_Neut_L_NAcc + .8*AWin_v_Neut_R_NAcc + .45*AWin_v_Neut_R_Insula +
            .7*BWin_v_Neut_L_NAcc +   .7*BWin_v_Neut_R_NAcc + .4*BWin_v_Neut_R_Insula +
            .8*BWin_v_BLose_L_NAcc +  .8*BWin_v_BLose_R_NAcc
                
Avoid =~  .8*ALose_v_Neut_L_Insula  + .8*ALose_v_Neut_R_Insula +
          .75*BLose_v_Neut_L_Insula + .75*BLose_v_Neut_R_Insula +
          .8*BLose_v_BWin_L_Insula  + .45*BLose_v_BWin_R_Insula
            
# Factor Covariances 
Approach ~~ -.6*Avoid

# Fixing factor variances
Approach ~~ 1*Approach
Avoid ~~ 1*Avoid

# Factor means/intercepts
#Approach ~ 1
#Avoid ~ 1
'

## General samples
### Using simsem to fit population model by creating a simulated set of 50 population sets. Do this for a dummy AHRB, MLS and ABCD sample
### the samples size is mean to be comparable to what I'd estimate I'd have access to in the real data.
set.seed(25151215)
sim_AHRB <- simsem::sim(nRep = 50, model = "lavaan", n = 104, 
                        generate = population_model, std.lv = TRUE, lavaanfun = "sem", 
                        # std.lv ~ ix the variances of all the latent variables 
                        dataOnly=T, meanstructure = FALSE, seed=123)

sim_MLS <- simsem::sim(nRep = 50, model = "lavaan", n = 120, 
                       generate = population_model, std.lv = TRUE, lavaanfun = "sem", 
                       dataOnly=T, meanstructure = FALSE, seed=123)

sim_ABCD <- simsem::sim(nRep = 50, model = "lavaan", n = 1000, 
                        generate = population_model, std.lv = TRUE, lavaanfun = "sem", 
                        dataOnly=T, meanstructure = FALSE, seed=123)


### for each simulate sets (50) of data, taking the mean of sets to create final study specific datasets, AHRB (3), MLS (2), ABCD (1)
sim_AHRB_data <- data.frame(aaply(laply(sim_AHRB, as.matrix), c(2,3), mean))
sim_AHRB_data$set <-3

sim_MLS_data <- data.frame(aaply(laply(sim_MLS, as.matrix), c(2,3), mean))
sim_MLS_data$set <-2

sim_ABCD_data <- data.frame(aaply(laply(sim_ABCD, as.matrix), c(2,3), mean))
sim_ABCD_data$set <-1

### Combined datasets
brain_set <- rbind(sim_AHRB_data,sim_MLS_data,sim_ABCD_data)

## Correlation matrix of data. 
###Using Hmisc to create a 24x24 matrix for a list (3) that contains: the pearson's r corr, sample size (N), and significance (p).
Brain_corr = rcorr(as.matrix(subset(brain_set,select=-c(set))), # excluding the set of data related to sample
                   type = "pearson")


### Using corrplot() to create heatmap of the data. 
par(mfrow=c(1,1))
corrplot(Brain_corr$r, type = "upper", 
         order = 'hclust',
         method =  "color", 
         tl.cex = 0.5, tl.col = 'black',
         cl.pos = 'r', tl.pos = 'lt', outline = TRUE,
         col=colorRampPalette(c("navyblue","white","red2"))(100),# colours http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
         mar = c(2,.15,.25,.15)#bottom, left, top and right,
)

# Running [restricted] CFA + Multigroup CFA 
## Specifiying apriori mdeol
MID_model <-'

# Factor loadings
Approach =~ AWin_v_Neut_L_NAcc_run1  + AWin_v_Neut_R_NAcc_run1  + AWin_v_Neut_R_Insula_run1 +
            BWin_v_Neut_L_NAcc_run1  + BWin_v_Neut_R_NAcc_run1  + BWin_v_Neut_R_Insula_run1 +
            BWin_v_BLose_L_NAcc_run1 + BWin_v_BLose_R_NAcc_run1 +
            AWin_v_Neut_L_NAcc_run2  + AWin_v_Neut_R_NAcc_run2 + AWin_v_Neut_R_Insula_run2 +
            BWin_v_Neut_L_NAcc_run2  + BWin_v_Neut_R_NAcc_run2  + BWin_v_Neut_R_Insula_run2 +
            BWin_v_BLose_L_NAcc_run2 + BWin_v_BLose_R_NAcc_run2 
                
Avoid =~    ALose_v_Neut_L_Insula_run1 + ALose_v_Neut_L_Insula_run1 +
            BLose_v_Neut_L_Insula_run1 + BLose_v_Neut_R_Insula_run1 +
            BLose_v_BWin_L_Insula_run1 + BLose_v_BWin_R_Insula_run1 +
            ALose_v_Neut_L_Insula_run2 + ALose_v_Neut_R_Insula_run2 +
            BLose_v_Neut_L_Insula_run2 + BLose_v_Neut_R_Insula_run2 +
            BLose_v_BWin_L_Insula_run2 + BLose_v_BWin_R_Insula_run2 
'



## Running CFA: Combined
all_sample <- cfa(model = MID_model, data = brain_set,
               estimator = "MLR", std.lv = TRUE, meanstructure = TRUE)

## Fitting Configural CFA
configural_cfa <- cfa(model = MID_model, data = brain_set, group = 'set', 
                      estimator = "MLR", std.lv = TRUE, meanstructure = TRUE)

## Fitting Metric CFA
metric_cfa <-cfa(model = MID_model, data = brain_set, 
                 group = 'set', group.equal=c("loadings"),
                 estimator = "MLR", std.lv = TRUE, meanstructure = TRUE)
## Extracting Fit Statistics
out <- matrix(NA, ncol = 9, nrow = 4)
colnames(out) <- c("model","chisq","df","pvalue", "rmsea", "cfi", "srmr",
                   "AIC", "BIC")

out[1,2:7] <- round(data.matrix(fitmeasures(all_sample, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi", "srmr"))), 
                    digits=3)

out[2,2:7] <- round(data.matrix(fitmeasures(configural_cfa, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi", "srmr"))), 
                            digits=3)

out[3,2:7] <- round(data.matrix(fitmeasures(metric_cfa, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi", "srmr"))), 
                    digits=3)


out[1,8] <- round(AIC(all_sample),3)
out[2,8] <- round(AIC(configural_cfa),3)
out[3,8] <- round(AIC(metric_cfa),3)
out[1,9] <- round(BIC(all_sample),3)
out[2,9] <- round(BIC(configural_cfa),3)
out[3,9] <- round(BIC(metric_cfa),3)

out[1:3,1] <-  c("Overall CFA", "Configg MG-CFA", "Metric MG-CFA")


## Model Parameter Summary
### All samples CFA
parameters(all_sample, standardize = T)

### Configural CFA model 
parameters(configural_cfa, standardize = T)

### Metric CFA model 
parameters(metric_cfa, standardize = T)

## Comparing models w/ BIC/AIC (anova)
anova(all_sample, configural_cfa)
anova(configural_cfa, metric_cfa)


# Running [semi-restricted] ESEM Model 
## ordering variables so can specify numerically
esem_data = brain_set[,c("AWin_v_Neut_L_NAcc_run1"  ,"AWin_v_Neut_L_NAcc_run2" ,
                         "BWin_v_Neut_L_NAcc_run1"  ,"BWin_v_Neut_L_NAcc_run2" ,
                         "BWin_v_BLose_L_NAcc_run1" ,"BWin_v_BLose_L_NAcc_run2",
                          "AWin_v_Neut_R_NAcc_run1" , "AWin_v_Neut_R_NAcc_run2",
                          "BWin_v_Neut_R_NAcc_run1" , "BWin_v_Neut_R_NAcc_run2",
                          "BWin_v_BLose_R_NAcc_run1", "BWin_v_BLose_R_NAcc_run2",
                         # insula values apprach 
                         "AWin_v_Neut_R_Insula_run1","AWin_v_Neut_R_Insula_run2", 
                         "BWin_v_Neut_R_Insula_run1","BWin_v_Neut_R_Insula_run2", 
                         # avoidance
                         "ALose_v_Neut_L_Insula_run1","ALose_v_Neut_L_Insula_run2",
                         "BLose_v_Neut_L_Insula_run1","BLose_v_Neut_L_Insula_run2",
                         "BLose_v_BWin_L_Insula_run1","BLose_v_BWin_L_Insula_run2",
                         "ALose_v_Neut_R_Insula_run1","ALose_v_Neut_R_Insula_run2",
                         "BLose_v_Neut_R_Insula_run1","BLose_v_Neut_R_Insula_run2",
                         "BLose_v_BWin_R_Insula_run1","BLose_v_BWin_R_Insula_run2",
                         "set")]
## Specify EFA Model; First, consistent w/ March et al. (2014), creating target rotation ensure they match onto variable list 
target_rot <- make_target(28,mainloadings = list(f1 = 1:16, f2 = 17:28))
esem.efa <- esem_efa(data = esem_data[,1:28], nfactors = 2,
                     target = target_rot, fm = "ml")
esem.efa$loadings

## per the example from Mateus Silverstrin, need to define anchor for each factor (value to loads highers on 1 factor and lowest on other)
anchor <- find_referents(efa_object = esem.efa,factor_names = c("f1","f2"))

## Pull starting parameters
esem_mid_model <- syntax_composer(efa_object = esem.efa, referents = anchor)


## Run ESEM model
### Specified Model

cat(esem_mid_model)

esem_mid_fit<- cfa(esem_mid_model, esem_data[,1:28], 
                   std.lv=TRUE, meanstructure = TRUE,
                   estimator = "MLR")

out[4,2:7] <- round(data.matrix(fitmeasures(esem_mid_fit, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi", "srmr"))), 
                    digits=3)
out[4,8] <- round(AIC(esem_mid_fit),3)
out[4,9] <- round(BIC(esem_mid_fit),3)
out[4,1] <-  c("Overall ESEM")

out <- as.data.frame(out)

out %>% 
  knitr::kable(
    col.names = c("Model", "Chi-sq", "DF", "p value", "RMSEA", "CFI", "SRMR", "AIC", "BIC"),
    caption = "Fit statistics from MG CFA and ESEM models",
    booktabs = TRUE
    )

# Running EFA [Unrestricted] model
## Rec. # Factors

### Using alt calc parallel analysis
paran(x = esem_data[,1:28],
      iterations = 1000, quietly = FALSE, centile = 95, 
      status = FALSE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      seed = 100)

plot(nScree(x=esem_data[,1:28],model="factors"))

### Using BIC values 

rec_factors <- matrix(NA, ncol = 2, nrow = 20)
colnames(rec_factors) <- c("Nfactors","BIC")

for (f in 1:20) {
  test_fac <- fa(r = esem_data[,1:28],  #raw data  
            nfactors = f, 
            rotate = "promax")
  rec_factors[f,1] <- f
  rec_factors[f,2] <-test_fac$BIC
}

bic_fact = as.data.frame(rec_factors)

lowest_bic <- which.min(bic_fact$BIC)

bic_fact %>% 
  ggplot(aes(x = Nfactors, y = BIC)) +
  geom_line(colour = 'black', linetype = 'dashed') +
  geom_vline(xintercept = bic_fact$Nfactors[lowest_bic], colour = 'red')+
  theme_minimal()

## Run EFA on sample specific
### EFA ABCD
abcd_efadata = subset(esem_data %>% filter(set==1))

abcd_efa <- factanal(x = abcd_efadata[,1:28],  #raw data  
              factors = 2, rotation = "promax" # oblique rotation allow for non-orthogonal structure
              )
heatmaply(round(abcd_efa$loadings[,1:2],2) %>% print(sort = T),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                 low = "blue", 
                 high = "darkred", 
                 space = "Lab",
                 midpoint = 0, 
                 limits = c(-1, 1)
               ),
               dendrogram = "none",
               xlab = "", ylab = "", 
               main = "",
               margins = c(60,100,40,20),
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               hide_colorbar = FALSE,
               branches_lwd = 0.1,
               label_names = c("Brain:", "Feature:", "Value"),
               fontsize_row = 9, fontsize_col = 9,
               labCol = colnames(abcd_efa$loadings[,1:2]),
               labRow = rownames(abcd_efa$loadings[,1:2]),
               heatmap_layers = theme(axis.line=element_blank()),
          
)

### EFA MLS

mls_efadata = subset(esem_data %>% filter(set==2))

mls_efa <- factanal(x = mls_efadata[,1:28],  #raw data  
              factors = 2, rotation = "promax" # oblique rotation allow for non-orthogonal structure
              )

heatmaply(round(mls_efa$loadings[,1:2],2) %>% print(sort = T),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                 low = "blue", 
                 high = "darkred", 
                 space = "Lab",
                 midpoint = 0, 
                 limits = c(-1, 1)
               ),
               dendrogram = "none",
               xlab = "", ylab = "", 
               main = "",
               margins = c(60,100,40,20),
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               hide_colorbar = FALSE,
               branches_lwd = 0.1,
               label_names = c("Brain:", "Feature:", "Value"),
               fontsize_row = 9, fontsize_col = 9,
               labCol = colnames(mls_efa$loadings[,1:2]),
               labRow = rownames(mls_efa$loadings[,1:2]),
               heatmap_layers = theme(axis.line=element_blank()),
          
)


### EFA AHRB

ahrb_efadata = subset(esem_data %>% filter(set==3))

ahrb_efa <- factanal(x = ahrb_efadata[,1:28],  #raw data  
              factors = 2, rotation = "promax" # oblique rotation allow for non-orthogonal structure
              )
heatmaply(round(ahrb_efa$loadings[,1:2],2) %>% print(sort = T),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
                 low = "blue", 
                 high = "darkred", 
                 space = "Lab",
                 midpoint = 0, 
                 limits = c(-1, 1)
               ),
               dendrogram = "none",
               xlab = "", ylab = "", 
               main = "",
               margins = c(60,100,40,20),
               grid_color = "white",
               grid_width = 0.00001,
               titleX = FALSE,
               hide_colorbar = FALSE,
               branches_lwd = 0.1,
               label_names = c("Brain:", "Feature:", "Value"),
               fontsize_row = 9, fontsize_col = 9,
               labCol = colnames(ahrb_efa$loadings[,1:2]),
               labRow = rownames(ahrb_efa$loadings[,1:2]),
               heatmap_layers = theme(axis.line=element_blank()),
          
)

## Comparing EFA models using congruence
fa.congruence(x = list(abcd_efa, mls_efa, ahrb_efa), digits = 2) %>% 
  knitr::kable(
    col.names = c("1. ABCD F1", "2. ABCD F2", "3. MLS F1", "4. MLS F2","5. AHRB F1", "6. AHRB F2"),
    caption = "ABCD, MLS and AHRB EFA Factor Congruence",
    booktabs = TRUE
  )

# Running Local SEM=
## Run LSEM
### first, creating a randomly generated pubertal scale variable to apply in simulated data
sim_ABCD_data$PDS <- as.integer(rtnorm(n=1000, mean = 3.5, sd = 1.5, 
                                       lower = 1, upper = 5))

# running LSEM
lsem.MID <- sirt::lsem.estimate(data = sim_ABCD_data, moderator = 'PDS', # moderator variable
                                moderator.grid = seq(1,5,1), # moderator levels, PDS 1 - 5
                                lavmodel = MID_model, # model
                                h = 2, # bandwidth parameter 
                                residualize = FALSE, # allow mean level differences 
                                meanstructure = TRUE,
                                std.lv=TRUE
)

## Summary LSEM
summary(lsem.MID)

## Plot LSEM
plot(lsem.MID, parindex=1:20)


## Permutation Test LSEM
lsem.permuted <- sirt::lsem.permutationTest(lsem.object = lsem.MID,
                                            B = 10, # permutations 
                                            residualize = FALSE) 
summary(lsem.permuted) # examine results