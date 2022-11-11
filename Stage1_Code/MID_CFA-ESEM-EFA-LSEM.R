# Automatically check and install necessary packages that were used in R v 4.0
if (!require("pacman")) install.packages("pacman")
pacman::p_load(simsem, plyr, tidyverse,ggdist, ggplot2, ggpubr, tidyquant, corrplot, reshape2,
               data.table,ggplot2,esemComp,
               Hmisc,nFactors,car,psych, paran, #for EFA
               semTools, semPlot, lavaan, # FOR CFA/ESEM
               parameters, superheat,weights, devtools, heatmaply, # For plotting/tabling
               sirt, msm
)


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

'


## General samples
set.seed(25151215)
sim_AHRB <- simsem::sim(nRep = 50, model = "lavaan", n = 108, 
                        generate = population_model, std.lv = TRUE, lavaanfun = "sem", 
                        # std.lv ~ ix the variances of all the latent variables 
                        dataOnly=T, meanstructure = FALSE, seed=123)

sim_MLS <- simsem::sim(nRep = 50, model = "lavaan", n = 159, 
                       generate = population_model, std.lv = TRUE, lavaanfun = "sem", 
                       dataOnly=T, meanstructure = FALSE, seed=123)

sim_ABCD <- simsem::sim(nRep = 50, model = "lavaan", n = 1000, 
                        generate = population_model, std.lv = TRUE, lavaanfun = "sem", 
                        dataOnly=T, meanstructure = FALSE, seed=123)

# sey pseudo site labels and combined
sim_AHRB_data <- data.frame(aaply(laply(sim_AHRB, as.matrix), c(2,3), mean))
sim_AHRB_data$set <-3

sim_MLS_data <- data.frame(aaply(laply(sim_MLS, as.matrix), c(2,3), mean))
sim_MLS_data$set <-2

sim_ABCD_data <- data.frame(aaply(laply(sim_ABCD, as.matrix), c(2,3), mean))
sim_ABCD_data$set <-1


brain_set <- rbind(sim_AHRB_data,sim_MLS_data,sim_ABCD_data)


## Correlation matrix of data
Brain_corr = rcorr(as.matrix(subset(brain_set,select=-c(set))), # excluding the set of data related to sample
                   type = "pearson")

par(mfrow=c(1,1))
corrplot(Brain_corr$r, type = "upper", 
         order = 'hclust',
         method =  "color", 
         tl.cex = 0.5, tl.col = 'black',
         cl.pos = 'r', tl.pos = 'lt', outline = TRUE,
         col=colorRampPalette(c("navyblue","white","red2"))(100),# colours http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
         mar = c(2,.15,.25,.15)#bottom, left, top and right,
)



# Running [restricted] Multigroup CFA 
  
## CFA model 
MID_model <-'

# Factor loadings
Approach =~ AWin_v_Neut_L_NAcc_run1  + AWin_v_Neut_R_NAcc_run1  + AWin_v_Neut_R_Insula_run1 +
            AWin_v_Neut_L_NAcc_run2  + AWin_v_Neut_R_NAcc_run2 + AWin_v_Neut_R_Insula_run2 +
            BWin_v_Neut_L_NAcc_run1  + BWin_v_Neut_R_NAcc_run1  + BWin_v_Neut_R_Insula_run1 +
            BWin_v_Neut_L_NAcc_run2  + BWin_v_Neut_R_NAcc_run2  + BWin_v_Neut_R_Insula_run2 +
            BWin_v_BLose_L_NAcc_run1 + BWin_v_BLose_R_NAcc_run1 +
            BWin_v_BLose_L_NAcc_run2 + BWin_v_BLose_R_NAcc_run2 
                
Avoid =~    ALose_v_Neut_L_Insula_run1 + ALose_v_Neut_R_Insula_run1 +
            ALose_v_Neut_L_Insula_run2 + ALose_v_Neut_R_Insula_run2 +
            BLose_v_Neut_L_Insula_run1 + BLose_v_Neut_R_Insula_run1 +
            BLose_v_Neut_L_Insula_run2 + BLose_v_Neut_R_Insula_run2 +
            BLose_v_BWin_L_Insula_run1 + BLose_v_BWin_R_Insula_run1 +
            BLose_v_BWin_L_Insula_run2 + BLose_v_BWin_R_Insula_run2 
'

## Running CFA: Three Samples
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
out <- matrix(NA, ncol = 10, nrow = 4)
colnames(out) <- c("model","chisq","df","pvalue", "rmsea", "cfi","tli", "srmr",
                   "AIC", "BIC")

out[1,2:8] <- round(data.matrix(fitmeasures(all_sample, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi","tli", "srmr"))), 
                    digits=3)

out[2,2:8] <- round(data.matrix(fitmeasures(configural_cfa, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi","tli", "srmr"))), 
                            digits=3)

out[3,2:8] <- round(data.matrix(fitmeasures(metric_cfa, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi","tli", "srmr"))), 
                    digits=3)


out[1,9] <- round(AIC(all_sample),3)
out[2,9] <- round(AIC(configural_cfa),3)
out[3,9] <- round(AIC(metric_cfa),3)

out[1,10] <- round(BIC(all_sample),3)
out[2,10] <- round(BIC(configural_cfa),3)
out[3,10] <- round(BIC(metric_cfa),3)

out[1:3,1] <-  c("Overall CFA", "Configg MG-CFA", "Metric MG-CFA")

## Model Parameter Summary


### All Sample CFA model 
parameters(all_sample, standardize = T)

### Configural CFA model 
parameters(configural_cfa, standardize = T)

### Metric CFA model 
parameters(metric_cfa, standardize = T)

## Comparing models w/ BIC/AIC (anova)
anova(all_sample, configural_cfa)

anova(configural_cfa, metric_cfa)
## Plotting multi-group config. CFA

layout(t(1:3))
semPaths(configural_cfa,
         color = "lightyellow",
         theme="colorblind",
         whatLabels = "std",
         style = "lisrel",
         sizeLat = 10,
         sizeLat2 = 10,
         sizeMan = 6,
         edge.color = "steelblue",
         edge.label.cex = 2,
         label.cex = 2,
         rotation = 2,
         layout = "tree2",
         intercepts = TRUE,
         residuals = FALSE,
         #residScale = 10,
         curve = 2,
         title = T,
         title.color = "black",
         cardinal = "lat cov",
         curvePivot = T,
         nCharNodes = 6,
         #nodeLabels = label,
         mar = c(2,5,2,6))
### Title 
title("Multi-group CFA on MID task Contrasts")



# Running [semi-restricted] ESEM Model {.tabset}

## Selected items for ESEM 

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



## Specify EFA Model 
target_rot <- make_target(28,mainloadings = list(f1 = 1:16, f2 = 17:28))
esem.efa <- esem_efa(data = esem_data[,1:28], nfactors = 2,
                     target = target_rot, fm = "ml")

esem.efa$loadings
anchor <- find_referents(efa_object = esem.efa,factor_names = c("f1","f2"))
esem_mid_model <- syntax_composer(efa_object = esem.efa, referents = anchor)

## Run ESEM model
### Specified Model

cat(esem_mid_model)

### Running full ESEM model 

esem_mid_fit<- cfa(esem_mid_model, esem_data[,1:28], std.lv=TRUE, meanstructure = TRUE,
                   estimator = "MLR")

out[4,2:8] <- round(data.matrix(fitmeasures(esem_mid_fit, 
                                            fit.measures = c("chisq","df","pvalue",
                                                             "rmsea", "cfi","tli", "srmr"))), 
                    digits=3)
out[4,9] <- round(AIC(esem_mid_fit),3)
out[4,10] <- round(BIC(esem_mid_fit),3)
out[4,1] <-  c("Overall ESEM")

out <- as.data.frame(out)

out %>% 
  knitr::kable(
    col.names = c("Model", "Chi-sq", "DF", "p value", "RMSEA", "CFI", "TLI","SRMR", "AIC", "BIC"),
    caption = "Fit statistics from MG-CFA and ESEM models",
    booktabs = TRUE
    )


# Aim: 2. Running EFA [Unrestricted] model 

## By Sample EFA

abcd_efa_df <- esem_data %>% filter(set == 1)
mls_efa_df <- esem_data %>% filter(set == 2)
ahrb_efa_df <- esem_data %>% filter(set == 3)


### Sample: *ABCD*

fa_abcd <- subset(abcd_efa_df[,1:28])

plot(nScree(x=fa_abcd,model="factors"))


paran(x = fa_abcd,
      iterations = 1000, quietly = FALSE, centile = 95, 
      status = FALSE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      seed = 100)


abcd_rec_factors <- matrix(NA, ncol = 2, nrow = 20)
colnames(abcd_rec_factors) <- c("Nfactors","BIC")

for (f in 1:20) {
  test_fac <- fa(r = fa_abcd,  #raw data  
            nfactors = f, fm = 'minres',
            rotate = "oblimin")
  abcd_rec_factors[f,1] <- f
  abcd_rec_factors[f,2] <-test_fac$BIC
}

abcd_bic_fact = as.data.frame(abcd_rec_factors)

abcd_lowest_bic <- which.min(abcd_bic_fact$BIC)

abcd_bic_fact %>% 
  ggplot(aes(x = Nfactors, y = BIC)) +
  geom_line(colour = 'black', linetype = 'dashed') +
  geom_vline(xintercept = abcd_bic_fact$Nfactors[abcd_lowest_bic], colour = 'red')+
  theme_minimal()



### Running efa
abcd_efa <- factanal(x = fa_abcd,  #raw data  
              factors = 2, fm ='minres', rotation = "promax" # oblique rotation allow for non-orthogonal structure
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

### Sample: *MLS*
fa_mls <- subset(mls_efa_df[,1:28])

plot(nScree(x=fa_mls,model="factors"))

paran(x = fa_mls,
      iterations = 1000, quietly = FALSE, centile = 95, 
      status = FALSE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      seed = 100)

mls_rec_factors <- matrix(NA, ncol = 2, nrow = 10)
colnames(mls_rec_factors) <- c("Nfactors","BIC")

for (f in 1:10) {
  test_fac <- fa(r = fa_mls,  #raw data  
            nfactors = f, fm = 'minres',
            rotate = "promax")
  mls_rec_factors[f,1] <- f
  mls_rec_factors[f,2] <-test_fac$BIC
}

mls_bic_fact = as.data.frame(mls_rec_factors)

mls_lowest_bic <- which.min(mls_bic_fact$BIC)

mls_bic_fact %>% 
  ggplot(aes(x = Nfactors, y = BIC)) +
  geom_line(colour = 'black', linetype = 'dashed') +
  geom_vline(xintercept = mls_bic_fact$Nfactors[mls_lowest_bic], colour = 'red')+
  theme_minimal()

### Running efa
mls_efa <- factanal(x = fa_mls,  #raw data  
              factors = 3, rotation = "promax" # oblique rotation allow for non-orthogonal structure
              )

heatmaply(round(mls_efa$loadings[,1:3],2) %>% print(sort = T),
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
               labCol = colnames(mls_efa$loadings[,1:3]),
               labRow = rownames(mls_efa$loadings[,1:3]),
               heatmap_layers = theme(axis.line=element_blank()),
          
)

### Sample: *AHRB*

fa_ahrb <- subset(ahrb_efa_df[,1:28])

plot(nScree(x=fa_ahrb,model="factors"))
paran(x = fa_ahrb,
      iterations = 1000, quietly = FALSE, centile = 95, 
      status = FALSE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      seed = 100)

ahrb_rec_factors <- matrix(NA, ncol = 2, nrow = 10)
colnames(ahrb_rec_factors) <- c("Nfactors","BIC")

for (f in 1:10) {
  test_fac <- fa(r = fa_ahrb,  #raw data  
            nfactors = f, fm = 'minres',
            rotate = "promax")
  ahrb_rec_factors[f,1] <- f
  ahrb_rec_factors[f,2] <-test_fac$BIC
}

ahrb_bic_fact = as.data.frame(ahrb_rec_factors)

ahrb_lowest_bic <- which.min(ahrb_bic_fact$BIC)

ahrb_bic_fact %>% 
  ggplot(aes(x = Nfactors, y = BIC)) +
  geom_line(colour = 'black', linetype = 'dashed') +
  geom_vline(xintercept = ahrb_bic_fact$Nfactors[ahrb_lowest_bic], colour = 'red')+
  theme_minimal()

ahrb_efa <- factanal(x = fa_ahrb,  #raw data  
              factors = 4, rotation = "promax" # oblique rotation allow for non-orthogonal structure
              )
### Running EFA
heatmaply(round(ahrb_efa$loadings[,1:4],2) %>% print(sort = T),
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
               labCol = colnames(ahrb_efa$loadings[,1:4]),
               labRow = rownames(ahrb_efa$loadings[,1:4]),
               heatmap_layers = theme(axis.line=element_blank()),
          
)

## Quant. Convergence

fa.congruence(x = list(abcd_efa, mls_efa, ahrb_efa), digits = 2) %>% 
  knitr::kable(
    col.names = c("1. ABCD F1", "2. ABCD F2", 
                  "3. MLS F1", "4. MLS F2","5. MLS F3",
                  "6. AHRB F1", "7. AHRB F2","8. AHRB F3", "8. AHRB F4"),
    caption = "ABCD, MLS and AHRB EFA Factor Congruence",
    booktabs = TRUE
  )

# Running Local SEM 

## Run LSEM
set.seed(110)
sim_ABCD_data$PDS <- as.integer(rtnorm(n=1000, mean = 3.5, sd = 1.5, 
                                       lower = 1, upper = 7))

lsem.MID <- sirt::lsem.estimate(data = sim_ABCD_data, moderator = 'PDS', # moderator variable
                                moderator.grid = seq(1,5,1), # moderator levels, PDS 1 - 5
                                lavmodel = MID_model, # model
                                h = 2, # bandwidth parameter 
                                residualize = FALSE, # allow mean level differences 
                                meanstructure = TRUE,
                                std.lv=TRUE
)


summary(lsem.MID)


## Plot LSEM
plot(lsem.MID, parindex=1:20)


## Permutation Test LSEM
lsem.permuted <- sirt::lsem.permutationTest(lsem.object = lsem.MID,
                                            B = 10, # permutations 
                                            residualize = FALSE) 
summary(lsem.permuted) # examine results

# Sensitivity Analyses 

## UM Specific CFA/ESEM/EFA

ABCD_site_spec <- sim_ABCD_data %>% filter(PDS == 3)

sim_AHRB_data$PDS <- NA
sim_MLS_data$PDS <- NA

UM_site_spec <- rbind(sim_AHRB_data,
                      sim_MLS_data,
                      ABCD_site_spec # UM specified ABCD data
)

### CFA

UM_all_sample <- cfa(model = MID_model, data = UM_site_spec,
                     estimator = "MLR", std.lv = TRUE, meanstructure = TRUE)

UM_config_cfa <- cfa(model = MID_model, data = UM_site_spec, group = 'set', 
                     estimator = "MLR", std.lv = TRUE, meanstructure = TRUE)

UM_metric_cfa <-cfa(model = MID_model, data = UM_site_spec, 
                    group = 'set', group.equal=c("loadings"),
                    estimator = "MLR", std.lv = TRUE, meanstructure = TRUE)

UMsite_out <- matrix(NA, ncol = 9, nrow = 3)
colnames(UMsite_out) <- c("model","chisq","df","pvalue", "rmsea", "cfi", "srmr",
                          "AIC", "BIC")


UMsite_out[1,2:7] <- round(data.matrix(fitmeasures(UM_all_sample, 
                                                   fit.measures = c("chisq","df","pvalue",
                                                                    "rmsea", "cfi", "srmr"))), 
                           digits=3)

UMsite_out[2,2:7] <- round(data.matrix(fitmeasures(UM_config_cfa, 
                                                   fit.measures = c("chisq","df","pvalue",
                                                                    "rmsea", "cfi", "srmr"))), 
                           digits=3)

UMsite_out[3,2:7] <- round(data.matrix(fitmeasures(UM_metric_cfa, 
                                                   fit.measures = c("chisq","df","pvalue",
                                                                    "rmsea", "cfi", "srmr"))), 
                           digits=3)

UMsite_out[1,8] <- round(AIC(UM_all_sample),3)
UMsite_out[2,8] <- round(AIC(UM_config_cfa),3)
UMsite_out[3,8] <- round(AIC(UM_metric_cfa),3)


UMsite_out[1,9] <- round(BIC(UM_all_sample),3)
UMsite_out[2,9] <- round(BIC(UM_config_cfa),3)
UMsite_out[3,9] <- round(BIC(UM_metric_cfa),3)

UMsite_out[1:3,1] <-  c("Overall CFA", "Config MG-CFA", "Metric MG-CFA")


UMsite_out %>% 
  knitr::kable(
    caption = "Fit statistics from MG-CFA and ESEM models",
    booktabs = TRUE
  )

### EFA

UM_abcd_efadata = subset(UM_site_spec[,c("AWin_v_Neut_L_NAcc_run1"  ,"AWin_v_Neut_L_NAcc_run2" ,
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
                                         "set")] %>% filter(set==1))
paran(x = UM_abcd_efadata[,1:28],
      iterations = 1000, quietly = FALSE, centile = 95, 
      status = FALSE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      seed = 100)
plot(nScree(x=UM_abcd_efadata[,1:28],model="factors"))

UM_abcd_efa <- factanal(x = UM_abcd_efadata[,1:28],  #raw data  
                        factors = 2, rotation = "promax" # oblique rotation allow for non-orthogonal structure
)

heatmaply(round(UM_abcd_efa$loadings[,1:2],2) %>% print(sort = T),
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
          labCol = colnames(UM_abcd_efa$loadings[,1:2]),
          labRow = rownames(UM_abcd_efa$loadings[,1:2]),
          heatmap_layers = theme(axis.line=element_blank()),
          
)


fa.congruence(x = list(UM_abcd_efa, mls_efa, ahrb_efa), digits = 2) %>% 
  knitr::kable(
    col.names = c("1. ABCD F1", "2. ABCD F2", 
                  "3. MLS F1", "4. MLS F2","5. MLS F3",
                  "6. AHRB F1", "7. AHRB F2","8. AHRB F3", "8. AHRB F4"),
    caption = "ABCD, MLS and AHRB EFA Factor Congruence",
    booktabs = TRUE
  )



## Resampling ABCD CFA/EFA

abcd_resamp_df <- subset(brain_set[,c("AWin_v_Neut_L_NAcc_run1"  ,"AWin_v_Neut_L_NAcc_run2" ,
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
                                      "set")] %>% filter(set==1))


set.seed(1111)
boot_cfa <- cfa(model = MID_model, data = abcd_resamp_df,
                estimator = "MLR", std.lv = TRUE, meanstructure = TRUE)

out_cfaboot <- bootstrapLavaan(object = boot_cfa, R = 1000, 
                               FUN = fitMeasures,
                               fit.measures=c("chisq","rmsea","cfi",
                                              "tli","srmr","AIC","BIC"),
                               parallel="multicore", ncpus=4)

cfaboot_df <- data.frame(out_cfaboot)


resampled_efa <- matrix(NA, ncol = 2, nrow = 1000) # resampling 50x for Config + Metric models
colnames(resampled_cfa) <- c("Sample","Factors")

samples = 1000

for (s in 1:samples) {
  
  sub_df <- sample_n(tbl = abcd_resamp_df, size = 1000, replace = TRUE)
  
  resampled_efa[s,1] <- s
  
  val <- nScree(x=sub_df[,1:28],model="factors")
  resampled_efa[s,2] <- as.integer(val$Components[3])
  
  
  
}

resampled_res <- data.frame(resampled_efa)

### plotting
n = 1000

ci_plt1 <- cfaboot_df %>% 
  gather(key = "FitIndex", value = "Statistic",
         srmr,rmsea) %>% 
  dplyr::group_by(FitIndex) %>% 
  dplyr::summarize(m = mean(Statistic), stdev = sd(Statistic)) %>% 
  ggplot(aes(x =FitIndex, y = m, fill = FitIndex, color=FitIndex)) +
  geom_point()+
  geom_errorbar(aes(ymin=m-(1.96*stdev/sqrt(n)),
                    ymax=m+(1.96*stdev/sqrt(n))), 
                width=.2,
                position=position_dodge(0.05))+
  labs(
    title = 'CFA: RMSEA & SRMR',
    x = 'Fit Stats',
    y = 'Type',
  )+
  theme_minimal()


ci_plt2 <- cfaboot_df %>% 
  gather(key = "FitIndex", value = "Statistic",
         cfi,tli) %>% 
  dplyr::group_by(FitIndex) %>% 
  dplyr::summarize(m = mean(Statistic), stdev = sd(Statistic)) %>% 
  ggplot(aes(x =FitIndex, y = m, fill = FitIndex, color=FitIndex)) +
  geom_point()+
  geom_errorbar(aes(ymin=m-(1.96*stdev/sqrt(n)),
                    ymax=m+(1.96*stdev/sqrt(n))), 
                width=.2,
                position=position_dodge(0.05))+
  labs(
    title = 'CFA: CFI & TLI',
    x = 'Fit Stats',
    y = 'Type',
  )+
  theme_minimal()


ci_plt3 <- cfaboot_df %>% 
  gather(key = "FitIndex", value = "Statistic",
         aic,bic) %>% 
  dplyr::group_by(FitIndex) %>% 
  dplyr::summarize(m = mean(Statistic), stdev = sd(Statistic)) %>% 
  ggplot(aes(x =FitIndex, y = m, fill = FitIndex, color=FitIndex)) +
  geom_point()+
  geom_errorbar(aes(ymin=m-(1.96*stdev/sqrt(n)),
                    ymax=m+(1.96*stdev/sqrt(n))), 
                width=.2,
                position=position_dodge(0.05))+
  labs(
    title = 'CFA: AIC & BIC',
    x = 'Fit Stats',
    y = 'Type',
  )+
  theme_minimal()


ci_plt4 = resampled_res %>%
  dplyr::summarize(m = mean(Factors), stdev = sd(Factors)) %>% 
  ggplot(aes(x ="", y = m)) +
  geom_point()+
  geom_errorbar(aes(ymin=m-(1.96*stdev/sqrt(n)),
                    ymax=m+(1.96*stdev/sqrt(n))), 
                width=.2,
                position=position_dodge(0.05))+
  ylim(0,5)+
  labs(
    title = 'EFA Factors',
    x = '',
    y = 'Number of Factors',
  )+
  theme_minimal()

ci_plt1;ci_plt2;ci_plt3;ci_plt4



### Plotting the distribution of values
plt1 <- cfaboot_df %>% 
  gather(key = "FitIndex", value = "Statistic",
         srmr,rmsea) %>% 
  ggplot(aes(x =FitIndex, y = Statistic, fill = FitIndex, color=FitIndex)) +
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, 
                       point_colour = NA, alpha = .5) +
  geom_boxplot(width = .2, outlier.shape = NA, alpha = .3) + 
  geom_jitter(width = .05, alpha = .5) +
  theme_minimal()+
  labs(
    title = 'CFA: RMSEA & SRMR',
    x = 'Fit Stats',
    y = 'Type',
  )

plt2 = cfaboot_df %>% 
  gather(key = "FitIndex", value = "Statistic",
         cfi,tli) %>% 
  ggplot(aes(x =FitIndex, y = Statistic, fill = FitIndex, color=FitIndex)) +
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, 
                       point_colour = NA, alpha = .5) +
  geom_boxplot(width = .2, outlier.shape = NA, alpha = .3) + 
  geom_jitter(width = .05, alpha = .5) +
  theme_minimal()+
  labs(
    title = 'CFA: CFI & TLI',
    x = 'Fit Stats',
    y = 'Type',
  )


plt3 = cfaboot_df %>% 
  gather(key = "FitIndex", value = "Statistic",
         aic,bic) %>% 
  ggplot(aes(x =FitIndex, y = Statistic, fill = FitIndex, color=FitIndex)) +
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, 
                       point_colour = NA, alpha = .5) +
  geom_boxplot(width = .2, outlier.shape = NA, alpha = .3) + 
  geom_jitter(width = .05, alpha = .5) +
  theme_minimal()+
  labs(
    title = 'CFA: AIC & BIC',
    x = 'Fit Stats',
    y = 'Type',
  )


avg = round(mean(resampled_res$Factors),1)
minimum = min(resampled_res$Factors)
maximum = max(resampled_res$Factors)

plt4 = resampled_res %>%
  ggplot(aes(x ="", y = Factors)) +
  ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, 
                       point_colour = NA, alpha = .5) +
  geom_jitter(width = .05, alpha = .5) +
  theme_minimal()+
  labs(
    title = 'EFA: Parallel Analysis Recommended Factors',
    subtitle = paste("Mean: ",avg," [Min: ", minimum, "Max:", maximum,"]"),
    x = 'Fit Stats',
    y = 'Type',
  )

plt1;plt2;plt3;plt4

