library("deSolve"); library("sensitivity"); library("parallel");library("deSolve"); library("parallel"); library("ggpubr"); library("reshape2")

rm(list=ls())
setwd("/Users/amorgan/Documents/PostDoc/Diff_Tax_Analysis/Theoretical_Analysis/Formalised_Analysis/Model_Output/Sensitivity_Anal")

# Import in the Data ------------------------------------------------------

FBD_import <- readRDS("Sens_Anal_FBD.RDS")
res_import <- readRDS("Sens_Anal_Res.RDS")

# Plotting the Output -----------------------------------------------------

#Infections

f99_dataframe_fbd <- data.frame(X=colnames(FBD_import$X), FBD_import$D1/FBD_import$V, 
                                1 - FBD_import$Dt/FBD_import$V - FBD_import$D1/FBD_import$V , 1 - FBD_import$Dt/FBD_import$V)
colnames(f99_dataframe_fbd)[-1] <- c("first.order","higher.order", "total.order")

plotdata_FBD <- melt(f99_dataframe_fbd, id.vars = "X", measure.vars = c("higher.order","first.order"))

plotdata_FBD$X <- factor(plotdata_FBD$X, levels = reorder(unique(plotdata_FBD$X), - f99_dataframe_fbd$total.order))

#Resistance

f99_dataframe_res <- data.frame(X=colnames(res_import$X), res_import$D1/res_import$V, 
                                1 - res_import$Dt/res_import$V - res_import$D1/res_import$V , 1 - res_import$Dt/res_import$V)
colnames(f99_dataframe_res)[-1] <- c("first.order","higher.order", "total.order")

plotdata_res <- melt(f99_dataframe_res, id.vars = "X", measure.vars = c("higher.order","first.order"))

plotdata_res$X <- factor(plotdata_res$X, levels = reorder(unique(plotdata_res$X), - f99_dataframe_res$total.order))


# Plotting the eFAST ------------------------------------------------------

#Infections
p_efast_fbd <- ggplot(plotdata_FBD, aes(fill = variable, x =  X, y = value)) + theme_bw() + 
  geom_col(color = "black",position= "stack") + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(legend.position=c(0.25, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(labels = c("Second Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = bquote(bold("eFAST total/first order sensitivity indices- Infections")),
       x ="", y = "Sensitivity Index")

#Resistance
p_efast_res <- ggplot(plotdata_res, aes(fill = variable, x =  X, y = value)) + theme_bw() + 
  geom_col(color = "black", position= "stack") + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(legend.position=c(0.25, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  scale_fill_manual(labels = c("Second Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = bquote(bold("eFAST total/first order sensitivity indices - Average Resistance")),
       x ="", y = "Sensitivity Index")

