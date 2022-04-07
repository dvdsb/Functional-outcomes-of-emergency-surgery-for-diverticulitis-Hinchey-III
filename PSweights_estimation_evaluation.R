gc()
options("install.lock"=FALSE)
options(scipen = 6, digits = 4)

library(tidyverse)
library(twang)
library(xtable)
library(see)
library(grid)
library(ggplotify)

path <- c("//home.gu.gu.se/home-XB$/xbocda/Documents/RCVR/dat/")
A <- readRDS(paste0(path, "laplav_data_201015.rds"))

A2 <- A %>% filter(Surgery != "Other") %>%
  mutate(Surgery = factor(Surgery, labels = c("Hartmann´s procedure" , "Laparoscopic lavage")))

A3 <- A2 %>% dplyr::select(studienummer , Surgery , age , kol , hjart_karl_sjd , immunosuppr ,  sepsis_op, CCI )

A4 <- A3 %>%
  mutate(kol = factor(case_when(kol == "No" ~"No/Don´t" , kol == "Yes" ~"Yes",
                                kol == "Don´t know" ~"No/Don´t")) ,
         hjart_karl_sjd = factor(case_when(hjart_karl_sjd == "No" ~"No/Don´t" , hjart_karl_sjd == "Yes" ~"Yes",
                                           hjart_karl_sjd == "Don´t know" ~"No/Don´t")) ,
         immunosuppr = factor(case_when(immunosuppr == "No" ~"No/Don´t" , immunosuppr == "Yes" ~"Yes",
                                        immunosuppr == "Don´t know" ~"No/Don´t")),
         sepsis_op = factor(case_when(sepsis_op == "No" ~"No/Don´t" , sepsis_op == "Yes" ~"Yes",
                                      sepsis_op == "Don´t know" ~"No/Don´t")),
         CCI_dik = case_when(CCI <= 32.5 ~0 , CCI > 32.5 ~1  ),
         age_cent = (age -  mean(age, na.rm = T))/sd(age, na.rm=T) )

A4 <- A4 %>% mutate(kol = as.numeric(kol)-1    ,
                    CVD = as.numeric(hjart_karl_sjd) -1 ,
                    immunosuppr  = as.numeric(immunosuppr )-1 ,
                    sepsis_op  = as.numeric( sepsis_op ) -1,
                    treatment = as.numeric(Surgery) -1  )
#___________________________________________________________________________________________________
# Estimate the propensity score
set.seed(7678463)

ps <-ps(   treatment ~
             age_cent  +
             kol +
             CVD +
             immunosuppr +
             sepsis_op ,
           data=as.data.frame(A4),
           interaction.depth = 3,
           shrinkage = 0.01,
           estimand="ATE",
           verbose=FALSE,
           stop.method=c("es.mean" ),
           n.trees=50000  )
#--------------------------------------------------------------------------------
# Evaluate the extent to which the score and it´s weighting achives balance 
# in the synthetic groups

names(ps)
ps$balance
bal.table(ps)

bal.tU <- as_tibble(bal.table(ps)$unw) %>% dplyr::select(tx.mn, ct.mn , std.eff.sz , p,  ks.pval)
bal.tW <- as_tibble(bal.table(ps)$es.mean.ATE) %>% dplyr::select(tx.mn, ct.mn , std.eff.sz , p,  ks.pval)

bal <- bal.table(ps)

pretty.tab <- bal$es.mean.ATE[,c("tx.mn","ct.mn","ks", "p" , "std.eff.sz")]
pretty.tab <- cbind(pretty.tab, bal$unw[,c("tx.mn","ct.mn","ks", "p" , "std.eff.sz")])
names(pretty.tab) <- c("E(Y1|t=1)","E(Y0|t=1)","KS","p" , "StDiff",
                       "E(Y1|t=1)","E(Y0|t=1)","KS","p" , "StDiff" )
bal2 <- xtable(pretty.tab,
               caption = "Balance of the treatment and comparison groups",
               label = "tab:balance",
               digits = c(0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2) ,
               align=c("l","r","r","r","r" , "r","r","r","r", "r", "r"))

bal2
#____________________________________________________________________________________
vrs <- c(1, 2, 3, 4, 5)
UNW <- as_tibble(bal.table(ps)$unw) %>%  dplyr::select(std.eff.sz ) %>%
  mutate(X = 1)
WN <- as_tibble(bal.table(ps)$es.mean.ATE) %>% dplyr::select(std.eff.sz ) %>%
  mutate(X = 2)
UNWW <- full_join(UNW, WN)
UNWW2 <- tibble(cbind(vrs, UNWW)) %>% mutate(abs = abs(std.eff.sz )) %>%
  mutate(X = factor(X, labels = c("Unweighted", "Weighted")),
         vars = factor(vrs, labels =  c("Age", "COPD", "Cardiovascular disease",
                                        "Immunosuppression", "Sepsis at surgery")))
ddd <-  UNWW2 %>%   group_by(X) %>%
  summarise( mean = mean(abs)) ; ddd
# Plot
PLT <-    ggplot(UNWW2, aes(X, abs , group = vars, lty = vars )) +
  geom_line(  size = 1) +
  theme_bw() +
  labs(   y = "Average absolute standardized effect size" ) +
  theme(legend.position = c(0.7, 0.8)) +
  theme (axis.title.y = element_text(size = 16)) +
  theme (axis.title.x = element_text(size = 12)) +
  theme(legend.title=element_blank()) +
  theme(legend.text = element_text( size = 14)) +
  theme(axis.text.x = element_text(face="bold",
                                   size=14, angle=0),
        axis.text.y = element_text(face="bold",
                                   size=14, angle=0)) +
  scale_x_discrete(name = "", limits = c("Unweighted", "Weighted"))

windows()
PLT
#________________________________________________________________________________________________________________________
# Plots
PS_fig_1 <- as.grob(plot(ps, plots=2))
PS_fig_2 <- as.grob(plot(ps, plots=3))
PS_fig_3 <- as.grob(plot(ps, plots=4))
PS_fig_4 <- as.grob(plot(ps, plots=5))

windows()
PS_fig_1234 <- see::plots(PS_fig_3, PS_fig_4  )
PS_fig_1234

#Whereas propensity score stratification requires considerable overlap in these spreads, excellent covariate balance can
#often be achieved with weights, even when the propensity scores estimated for the treatment and
#control groups show little overlap, that is lack of overlap is compensated by weighting."

SMM0 <- summary(ps)
SMM0
#Variability decrease by weighting. Effective Sample Size (ESS) demonstrated how many "effective" patients
#we have left after weighting. It decrease from  173 to 141 and from 291 to 265 for treatment and control, respectively.
#This is the price we pay in terms of loss of precision in order to handle the bias due to confounding.
SMM1 <- summary(ps$gbm.obj,  plot=F)
SMM1
# Relative influence på propensity score: Age and immunosuppressive therapy had createst infleunce (56% resp 22%).

w  <- get.weights(ps,  stop.method = "es.mean")

anadata.w <- cbind(A4,w)
dim(anadata.w)

# Save data
write_rds(anadata.w , "S://RCVR/kod/LapLavQuestionnaire/LapLav_PS_w_ATE.rds")

