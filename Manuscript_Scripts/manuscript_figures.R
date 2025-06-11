library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ggsci)
library(viridis)
library(pbapply)
library(GWASTools)
library(gridExtra)
library(xtable)
library(readr)
library(VennDiagram)
library(writexl)
library(UpSetR)
library(tidyr)


###########################################################################
#### child functions ######################################################
read_in_r2 <- function(fname,dir){
  print(fname)
  tmp <- fread(paste0(dir,"/results/",fname),header=F,sep="\t",fill=T)
  tmp <- tmp[,-1]
  colnames(tmp) <- c("R2","N","Model","Sim","Prop_Causal_Grab","N_Causal_Grab_AFR","N_Causal_Grab_EUR")
  param_row <- as.numeric(strsplit(strsplit(fname,"_")[[1]][2],"row")[[1]][2])
  tmp$param_row <- param_row
  return(tmp)
}

read_in_twas <- function(fname,dir){
  print(fname)
  tmp <- read.table(paste0(dir,fname),header=T,sep="\t")
  return(tmp)
}

read_in_snp_cnt <- function(fname,dir){
  print(fname)
  tmp <- read.table(paste0(dir,"/results/snp_counts/",fname),header=T,sep="\t")
  pr <- as.numeric(strsplit(fname,"_")[[1]][1])
  tmp$pr <- pr
  return(tmp)
}
###########################################################################
###########################################################################

pass <- 36
dir <- paste0("/Volumes/walden/Admix_TWAS_Project/pass",pass)
# set working directory for figure output
out_dir="/Volumes/walden/Admix_TWAS_Project/manuscript/v8/figures/"

# read in parameter space
parameter_space <- read.table(paste0(dir,"/files_for_analysis/paramspace_full.txt"),sep="\t",header=T)
parameter_space$param_row <- 1:nrow(parameter_space)

# read in result file names
file.list <- list.files(paste0(dir,"/results"),recursive = T)
twas.files <- file.list[grep("acat_res",file.list)]
r2.files <- file.list[grep("R2",file.list)]
h2.files <- file.list[grep("he2",file.list)]

# process imputation R2 files
out <- lapply(r2.files,FUN=read_in_r2,dir=dir)
dat <- bind_rows(out, .id = "column_label")

# clean up dataframe structure
dat <- merge(dat,parameter_space,by="param_row")
dat$PRS <- unlist(lapply(strsplit(dat$Model,"_"), `[[`, 1))
dat$PRS <- factor(dat$PRS)
dat$rho <- factor(dat$rho,levels=c(0.5,1))
dat$causal_num <- factor(dat$causal_num,levels=unique(dat$causal_num))
dat$h2 <- paste0("H2 AFR ",dat$he2_afr, "/H2 EUR ", dat$he2_eur)
dat$h2 <- factor(dat$h2)
dat$Testing_N <- factor(dat$N)

# exclude weighted casPS simulation results
dat <- dat[-grep("casPS_wt",dat$Model),]
dat$Method <- NA
dat$Method[grep("casPS",dat$Model)] <- "casPS"
dat$Method[grep("AFR_500_aspPS",dat$Model)] <- "aspPS AFR"
dat$Method[grep("EUR_500_aspPS",dat$Model)] <- "aspPS EUR"
dat$Method[grep("ADMIX_80_10g_500",dat$Model)] <- "ADMIX_80_10g PS"
dat$Method[grep("ADMIX_80_5g_500",dat$Model)] <- "ADMIX_80_5g PS"
dat$Method[grep("ADMIX_50_5g_500",dat$Model)] <- "ADMIX_50_5g PS"
dat$Method[grep("ref_EUR_500_standard PS",dat$Model)] <- "EUR PS"
dat$Method[grep("ref_AFR_500_standard PS",dat$Model)] <- "AFR PS"
dat$Method <- factor(dat$Method,levels=c("aspPS AFR",
                                         "aspPS EUR",
                                         "casPS",
                                         "ADMIX_80_10g PS",
                                         "ADMIX_80_5g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"))
dat$op <- factor(dat$op,levels=unique(dat$op))
dat$Misclass_Pct <- "0%"
dat$Misclass_Pct[grep("mcpt_10",dat$Model)] <- "10%"
dat$Misclass_Pct <- factor(dat$Misclass_Pct)
dat$Approach <- "Standard Imputation"
dat$Approach[dat$Method %in% c("aspPS AFR","aspPS EUR","casPS")] <- "LA-Aware Imputation"

model_vals <- unique(dat$Model)

dat$group_full <- paste(dat$h2,
                        dat$rho,
                        dat$Testing_N,
                        dat$Misclass_Pct,
                        dat$op,
                        dat$causal_num,
                        dat$Method,
                        dat$PRS)
whisker <- tapply(dat$R2,INDEX=dat$group_full,FUN=function(x){
  iqr <- IQR(x)
  q3 <- quantile(x)[4]
  x <- x[x<=(1.5*iqr+q3)]
  return(max(x))
})
max_y <- max(whisker)

op_names <- list(
  '0.5'="OP 0.5",
  '1'="OP 1"
)
causal_num_names <- list(
  '2'="2 eQTLs",
  '10'="10 eQTLs",
  '100'="100 eQTLs"
)

testing_n_names <- list(
  '1000'="Test N 1000",
  '5000'="Test N 5000",
  '10000'="Test N 10000"
)
rho_names <- list(
  '0.5'= "Rho 0.5",
  '1'= "Rho 1"
)

plot_labeller <- function(variable,value){
  if (variable=='op') {
    return(op_names[value])
  } else if (variable=='causal_num') {
    return(causal_num_names[value])
  } else if (variable=='Testing_N') {
    return(testing_n_names[value])
  } else if (variable=='rho') {
    return(rho_names[value])
  } else {
    return(as.character(value))
  }
}



# imputation R2
# 2 eQTL
# He2 afr 0.2 / eur 0.2
# no LA misclass


tmp <- dat %>%
  filter(causal_num==2, he2_eur == 0.2, he2_afr==0.2, Misclass_Pct=="0%")

tmp <- tmp[!is.na(tmp$Method),]
tmp$rho_label <- paste0("rho == ", tmp$rho)
tmp$op <- paste0("OP = ",tmp$op)
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=PRS)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2")+theme_bw()+
  ylab(expression(paste("Imputation ",R^{2})))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill="PRS Method",title=expression("2 eQTLs, "~h[e]^2~AFR~"= 0.2,"~h[e]^2~EUR~"= 0.2"),
       subtitle = "LA Misclassification 0%")+
  facet_grid(rho_label~op, labeller = labeller(rho_label = label_parsed))+
  scale_y_continuous(limits=c(0,0.3))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(face=c(rep("bold",3),rep("plain",4)),angle = 90, vjust = 0.5, hjust=1, color = "black"),
        legend.text=element_text(size = 13),
        legend.title=element_text(size = 13))+scale_fill_manual(values=c("gold",
                                                                         "Deep Sky Blue 3",
                                                                         "Medium Sea Green"))
pdf(file=paste0(out_dir,"imp_r2_2eqtl_he2_afr_0.2_eur_0.2_no_LA_misclass.pdf"),
    width=10, height=8)
print(gg)
dev.off()


# imputation R2
# 10 eQTL
# He2 afr 0.2 / eur 0.2
# no LA misclass

tmp <- dat %>%
  filter(causal_num==10, he2_eur == 0.2, he2_afr==0.2, Misclass_Pct=="0%")

tmp <- tmp[!is.na(tmp$Method),]
tmp$rho_label <- paste0("rho == ", tmp$rho)
tmp$op <- paste0("OP = ",tmp$op)
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=PRS)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2")+theme_bw()+
  ylab(expression(paste("Imputation ",R^{2})))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill="PRS Method",title=expression("10 eQTLs, "~h[e]^2~AFR~"= 0.2,"~h[e]^2~EUR~"= 0.2"),
       subtitle = "LA Misclassification 0%")+
  facet_grid(rho_label~op, labeller = labeller(rho_label = label_parsed))+
  scale_y_continuous(limits=c(0,0.3))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(face=c(rep("bold",3),rep("plain",4)),angle = 90, vjust = 0.5, hjust=1, color = "black"),
        legend.text=element_text(size = 13),
        legend.title=element_text(size = 13))+scale_fill_manual(values=c("gold",
                                                                         "Deep Sky Blue 3",
                                                                         "Medium Sea Green"))
pdf(file=paste0(out_dir,"imp_r2_10eqtl_he2_afr_0.2_eur_0.2_no_LA_misclass.pdf"),
    width=10, height=8)
print(gg)
dev.off()


# supplemental
# imputation R2
# all eQTL
# He2 afr 0.2 / eur 0.2
# no LA misclass

tmp <- dat %>%
  filter(he2_eur == 0.2, he2_afr==0.2, Misclass_Pct=="0%")

tmp <- tmp[!is.na(tmp$Method),]
tmp$rho_label <- paste0("rho == ", tmp$rho)
tmp$op <- paste0("OP = ",tmp$op)
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10,100),labels=paste0(c(2,10,100)," eQTLs"))
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=PRS)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2")+theme_bw()+
  ylab(expression(paste("Imputation ",R^{2})))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill="PRS Method",title=expression(~h[e]^2~AFR~"= 0.2,"~h[e]^2~EUR~"= 0.2"),
       subtitle = "LA Misclassification 0%")+
  facet_grid(op+rho_label~causal_num, labeller = labeller(rho_label = label_parsed))+
  scale_y_continuous(limits=c(0,0.3))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(face=c(rep("bold",3),rep("plain",4)),angle = 90, vjust = 0.5, hjust=1, color = "black"),
        legend.text=element_text(size = 13),
        legend.title=element_text(size = 13))+scale_fill_manual(values=c("gold",
                                                                         "Deep Sky Blue 3",
                                                                         "Medium Sea Green"))
pdf(file=paste0(out_dir,"sfig_imp_r2_h2_afr_0.2_eur_0.2_no_LA_misclass.pdf"),
    width=8, height=10)
print(gg)
dev.off()


# supplemental
# imputation R2
# all eQTL
# He2 afr 0.2 /eur 0.1
# no LA misclass

tmp <- dat %>%
  filter(he2_eur == 0.1, he2_afr==0.2, Misclass_Pct=="0%")

tmp <- tmp[!is.na(tmp$Method),]
tmp$rho_label <- paste0("rho == ", tmp$rho)
tmp$op <- paste0("OP = ",tmp$op)
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10,100),labels=paste0(c(2,10,100)," eQTLs"))
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=PRS)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2")+theme_bw()+
  ylab(expression(paste("Imputation ",R^{2})))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill="PRS Method",title=expression(~h[e]^2~AFR~"= 0.2,"~h[e]^2~EUR~"= 0.1"),
       subtitle = "LA Misclassification 0%")+
  facet_grid(op+rho_label~causal_num, labeller = labeller(rho_label = label_parsed))+
  scale_y_continuous(limits=c(0,0.3))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(face=c(rep("bold",3),rep("plain",4)),angle = 90, vjust = 0.5, hjust=1, color = "black"),
        legend.text=element_text(size = 13),
        legend.title=element_text(size = 13))+scale_fill_manual(values=c("gold",
                                                                         "Deep Sky Blue 3",
                                                                         "Medium Sea Green"))
pdf(file=paste0(out_dir,"sfig_imp_r2_h2_afr_0.2_eur_0.1_no_LA_misclass.pdf"),
    width=8, height=10)
print(gg)
dev.off()

# supplemental
# imputation R2
# all eQTL
# He2 afr 0.1 / eur 0.2
# no LA misclass


tmp <- dat %>%
  filter(he2_eur == 0.2, he2_afr==0.1, Misclass_Pct=="0%")

tmp <- tmp[!is.na(tmp$Method),]
tmp$rho_label <- paste0("rho == ", tmp$rho)
tmp$op <- paste0("OP = ",tmp$op)
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10,100),labels=paste0(c(2,10,100)," eQTLs"))
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=PRS)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2")+theme_bw()+
  ylab(expression(paste("Imputation ",R^{2})))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(fill="PRS Method",title=expression(~h[e]^2~AFR~"= 0.1,"~h[e]^2~EUR~"= 0.2"),
       subtitle = "LA Misclassification 0%")+
  facet_grid(op+rho_label~causal_num, labeller = labeller(rho_label = label_parsed))+
  scale_y_continuous(limits=c(0,0.3))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(face=c(rep("bold",3),rep("plain",4)),angle = 90, vjust = 0.5, hjust=1, color = "black"),
        legend.text=element_text(size = 13),
        legend.title=element_text(size = 13))+scale_fill_manual(values=c("gold",
                                                                         "Deep Sky Blue 3",
                                                                         "Medium Sea Green"))
pdf(file=paste0(out_dir,"sfig_imp_r2_h2_afr_0.1_eur_0.2_no_LA_misclass.pdf"),
    width=8, height=10)
print(gg)
dev.off()

# supplemental
# imputation R2 with misclass
# P0.001


cols <- viridis(6)
cols_1 <- cols[c(1,3,5)]

prs="P0.001"
tmp <- dat[dat$Method%in%c("casPS","aspPS AFR","aspPS EUR"),]
tmp <- tmp[tmp$Testing_N=="10000",]
tmp <- tmp[!tmp$causal_num=="100",]
tmp <- tmp[tmp$PRS==prs,]
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=h2,alpha=Misclass_Pct)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2",title=paste0(prs))+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_manual(values=cols_1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_alpha_manual(values=c(1,0.5))+
  scale_color_manual(values=c("black","gold"))+
  #scale_color_manual(values=c("black","gold"))+
  labs(fill="Heritability",alpha=NULL)+theme_bw()+
  facet_grid(rho+op~causal_num,labeller=plot_labeller)+scale_y_continuous(limits=c(0,0.26))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=13),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(alpha = FALSE)
pdf(file=paste0(out_dir,"sfig_imp_R2_with_misclass_P0.001.pdf"),
    width=8, height=8)
print(gg)
dev.off()

# supplemental
# imputation R2 with misclass
# P0.05

prs="P0.05"
tmp <- dat[dat$Method%in%c("casPS","aspPS AFR","aspPS EUR"),]
tmp <- tmp[tmp$Testing_N=="10000",]
tmp <- tmp[!tmp$causal_num=="100",]
tmp <- tmp[tmp$PRS==prs,]
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=h2,alpha=Misclass_Pct)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2",title=paste0(prs))+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_manual(values=cols_1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_alpha_manual(values=c(1,0.5))+
  scale_color_manual(values=c("black","gold"))+
  #scale_color_manual(values=c("black","gold"))+
  labs(fill="Heritability",alpha=NULL)+theme_bw()+
  facet_grid(rho+op~causal_num,labeller=plot_labeller)+scale_y_continuous(limits=c(0,0.26))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=13),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(alpha = FALSE)
pdf(file=paste0(out_dir,"sfig_imp_R2_with_misclass_P0.05.pdf"),
    width=8, height=8)
print(gg)
dev.off()

# supplemental
# imputation R2 with misclass
# lassosum

prs="lassosum"
tmp <- dat[dat$Method%in%c("casPS","aspPS AFR","aspPS EUR"),]
tmp <- tmp[tmp$Testing_N=="10000",]
tmp <- tmp[!tmp$causal_num=="100",]
tmp <- tmp[tmp$PRS==prs,]
gg <- ggplot(tmp, aes(x=Method, y=R2, fill=h2,alpha=Misclass_Pct)) +geom_boxplot(outlier.shape = NA)+
  labs(x="",y="Testing R2",title=paste0(prs))+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #scale_fill_manual(values=cols_1)+
  scale_fill_brewer(palette = "Dark2")+
  scale_alpha_manual(values=c(1,0.5))+
  scale_color_manual(values=c("black","gold"))+
  #scale_color_manual(values=c("black","gold"))+
  labs(fill="Heritability",alpha=NULL)+theme_bw()+
  facet_grid(rho+op~causal_num,labeller=plot_labeller)+scale_y_continuous(limits=c(0,0.26))+
  theme(strip.text.x = element_text(size = 13),
        strip.text.y = element_text(size = 13),
        axis.text = element_text(size = 13),
        legend.text=element_text(size=11),
        legend.title=element_text(size=13),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(alpha = FALSE)
pdf(file=paste0(out_dir,"sfig_imp_R2_with_misclass_lassosum.pdf"),
    width=8, height=8)
print(gg)
dev.off()
  
####################################################################################
# qq plots
####################################################################################
dir <- "/Volumes/walden/Admix_TWAS_Project/pass36/results/"
out <- lapply(twas.files,FUN=read_in_twas,dir=dir)
dat <- bind_rows(out)
null <- dat[dat$hp2==0,]
alt <- dat[dat$hp2>0,]
rm(dat)

dir="/Volumes/walden/Admix_TWAS_Project/manuscript/v8/figures/"

# null acatl level 1, misclass 0, test N 10,000 AGGREGATE

tmp <- null[null$acat_level=="1: ACAT(prs)",]
tmp <- tmp[is.na(tmp$mcpt) | tmp$mcpt=="0%",]
tmp <- tmp[tmp$testing_n==10000,]
tmp$Method <- factor(tmp$m,levels=c("aspPS AFR",
                                    "aspPS EUR",
                                    "casPS",
                                    "ADMIX_80_10g PS",
                                    "ADMIX_80_5g PS",
                                    "ADMIX_50_5g PS",
                                    "AFR PS",
                                    "EUR PS"))

pdf(file=paste0(dir,"qq_acat_lev1_all_paramrow_N10000_mcpt0.pdf"),
      width=6, height=8)
par(mfrow=c(3,2))
for(method in c("aspPS AFR",
                  "aspPS EUR",
                  "casPS",
                  "ADMIX_80_10g PS",
                  "AFR PS",
                  "EUR PS")){
  tt <- tmp[tmp$Method==method,]
  qqPlot(tt$p,main=method)
}
dev.off()
par(mfrow=c(1,1))

# null acat level 1, misclass 10, test N 10,000 AGGREGATE

tmp <- null[null$acat_level=="1: ACAT(prs)",]
tmp <- tmp[tmp$mcpt=="10%",]
tmp <- tmp[tmp$testing_n==10000,]
tmp$Method <- factor(tmp$m,levels=c("aspPS AFR",
                                    "aspPS EUR",
                                    "casPS",
                                    "ADMIX_80_10g PS",
                                    "ADMIX_80_5g PS",
                                    "ADMIX_50_5g PS",
                                    "AFR PS",
                                    "EUR PS"))
pdf(file=paste0(dir,"qq_acat_lev1_all_paramrow_N10000_mcpt10.pdf"),
    width=8, height=4)
par(mfrow=c(1,3))
for(method in c("aspPS AFR",
                "aspPS EUR",
                "casPS")){
  tt <- tmp[tmp$Method==method,]
  qqPlot(tt$p,main=method)
}
dev.off()
par(mfrow=c(1,1))

# null acatl level 2, misclass 0, test N 10,000 AGGREGATE

tmp <- null[null$acat_level=="2: ACAT(ACAT(prs))",]
tmp <- tmp[tmp$testing_n==10000,]
tmp$m <- factor(tmp$m,levels=c("casPS,asPS",
                               "casPS,asPS,ref PS",
                               "casPS,asPS,ref PS,ADMIX_50_5g PS",
                               "casPS,asPS,ref PS,ADMIX_80_10g PS",
                               "casPS,asPS,ref PS,ADMIX_80_5g PS"),
                labels=c("casPS, aspPS",
                         "casPS, aspPS, ref PS",
                         "casPS, aspPS, ref PS, ADMIX_50_5g PS",
                         "casPS, aspPS, ref PS, ADMIX_80_10g PS",
                         "casPS, aspPS, ref PS, ADMIX_80_5g PS"))

pdf(file=paste0(dir,"qq_acat_lev2_all_paramrow_N10000.pdf"),
    width=8.5, height=6)
par(mfrow=c(2,3))

tmp2 <- tmp[tmp$mcpt=="0%",]
for(method in c("casPS, aspPS",
                "casPS, aspPS, ref PS",
                "casPS, aspPS, ref PS, ADMIX_80_10g PS")){
  tt <- tmp2[tmp2$m==method,]
  qqPlot(tt$p,main=paste0(method),sub="No LA Misclass",cex.main=0.93,cex.axis=1.2)
}
tmp2 <- tmp[tmp$mcpt=="10%",]
for(method in c("casPS, aspPS",
                "casPS, aspPS, ref PS",
                "casPS, aspPS, ref PS, ADMIX_80_10g PS")){
  tt <- tmp2[tmp2$m==method,]
  qqPlot(tt$p,main=paste0(method),sub="LA Misclass 10%",cex.main=0.93,cex.axis=1.2)
}
dev.off()
par(mfrow=c(1,1))

####################################################################################
# power figures
####################################################################################

# POWER alpha 2.5e-6
calc_power <- function(p,alpha){
  pow <- sum(p<alpha)/sum(!is.na(p))
  return(pow)
}
alpha <- 2.5e-6

parameter_space$h2 <- paste0("H2 AFR ",parameter_space$he2_afr, "/H2 EUR ", parameter_space$he2_eur)

tmp <- alt
tmp$id <- paste(tmp$param_row, tmp$testing_n, tmp$mcpt, tmp$m,sep="-")
t <- data.frame(tapply(tmp$p,INDEX=list(tmp$id),FUN=calc_power,alpha=alpha))
colnames(t) <- "power"
t$id <- row.names(t)
t$param_row <- as.numeric(unlist(lapply(strsplit(t$id,"-"),`[[`,1)))
t$testing_n <- as.numeric(unlist(lapply(strsplit(t$id,"-"),`[[`,2)))
t$mcpt <- unlist(lapply(strsplit(t$id,"-"),`[[`,3))
t$m <- unlist(lapply(strsplit(t$id,"-"),`[[`,4))
pdat <- merge(t,parameter_space,by="param_row")
pdat$method <- factor(pdat$m,levels=c("casPS,asPS",
                                      "casPS,asPS,ref PS",
                                      "casPS,asPS,ref PS,ADMIX_50_5g PS",
                                      'casPS,asPS,ref PS,ADMIX_80_10g PS',
                                      "casPS,asPS,ref PS,ADMIX_80_5g PS",
                                      "aspPS AFR",
                                      "aspPS EUR",
                                      "casPS",
                                      "ADMIX_80_10g PS",
                                      "ADMIX_80_5g PS",
                                      "ADMIX_50_5g PS",
                                      "AFR PS",
                                      "EUR PS"
),labels=c("casPS, aspPS",
           "casPS, aspPS, ref PS",
           "casPS, aspPS, ref PS, ADMIX_50_5g PS",
           'casPS, aspPS, ref PS, ADMIX_80_10g PS',
           "casPS, aspPS, ref PS, ADMIX_80_5g PS",
           "aspPS AFR",
           "aspPS EUR",
           "casPS",
           "ADMIX_80_10g PS",
           "ADMIX_80_5g PS",
           "ADMIX_50_5g PS",
           "AFR PS",
           "EUR PS"
))
pdat$Approach <- "1: ACAT(prs)"
pdat$Approach[pdat$method %in% c("casPS, aspPS",
                                 "casPS, aspPS, ref PS",
                                 "casPS, aspPS, ref PS, ADMIX_50_5g PS",
                                 'casPS, aspPS, ref PS, ADMIX_80_10g PS',
                                 "casPS, aspPS, ref PS, ADMIX_80_5g PS"
)] <- "2: ACAT(ACAT(prs))"
pdat$Approach <- factor(pdat$Approach)
pdat$LA_aware <- F
pdat$LA_aware[pdat$method %in% c("casPS, aspPS",
                                 "casPS, aspPS, ref PS",
                                 "casPS, aspPS, ref PS, ADMIX_50_5g PS",
                                 'casPS, aspPS, ref PS, ADMIX_80_10g PS',
                                 "casPS, aspPS, ref PS, ADMIX_80_5g PS",
                                 "casPS","aspPS AFR","aspPS EUR"
)] <- T
pdat$LA_aware <- factor(pdat$LA_aware)

# no misclass 
# 2 eqtls

tmp <- pdat[pdat$mcpt %in% c("0%","NA"),]
tmp <- tmp[tmp$testing_n==10000,]
tmp <- tmp[tmp$causal_num==2,]
tmp$rho_label <- paste0("rho == ", tmp$rho)
tmp$he2 <- paste0("afr",tmp$he2_afr,"eur",tmp$he2_eur)
tmp$he2 <- factor(tmp$he2,levels=c('afr0.1eur0.2', 'afr0.2eur0.1', 'afr0.2eur0.2'),
                  labels=c('He2 0.1/0.2',
                           'He2 0.2/0.1',
                           'He2 0.2/0.2'))
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10),labels=c("2 eQTLs","10 eQTLs"))
tmp$op <- factor(tmp$op,levels=c(0.5,1),labels=c("OP 0.5","OP 1"))
tmp <- tmp[tmp$method %in% c("casPS, aspPS, ref PS",
                             "aspPS AFR",
                             "aspPS EUR",
                             "casPS",
                             "ADMIX_80_10g PS",
                             "ADMIX_50_5g PS",
                             "AFR PS",
                             "EUR PS"),]
tmp$plot_method <- "Level 1 (LA-Unaware)"
tmp$plot_method[tmp$method %in% c("casPS, aspPS",
                                  "casPS, aspPS, ref PS")] <- "Level 2 (CADET)"
tmp$plot_method[tmp$method %in% c("aspPS AFR",
                                  "aspPS EUR",
                                  "casPS")] <- "Level 1 (LA-Aware)"
tmp$plot_method <- factor(tmp$plot_method,levels=c("Level 1 (LA-Unaware)",
                                                   "Level 1 (LA-Aware)",
                                                   "Level 2 (CADET)"))
tmp$method <- factor(tmp$method,levels=c("casPS, aspPS, ref PS",
                                         "aspPS AFR",
                                         "aspPS EUR",
                                         "casPS",
                                         "ADMIX_80_10g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"),
                     labels=c("CADET",
                              "aspPS AFR",
                              "aspPS EUR",
                              "casPS",
                              "ADMIX_80_10g PS",
                              "ADMIX_50_5g PS",
                              "AFR PS",
                              "EUR PS"))
tmp$rho_label <- paste0("rho == ", tmp$rho)
gg <- ggplot(data=tmp, aes(x=method, y=power,fill=plot_method)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
  labs(y=paste0("Power at alpha ",alpha),title="2 eQTLs, No LA Misclassification",subtitle="Testing N 10000",fill="Method",x="")+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black"),
        legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        legend.position="bottom")+
  coord_flip()+
  facet_grid(he2~op+rho_label, labeller = labeller(rho_label = label_parsed))+
  scale_fill_manual(values=c("Orchid","Deep Sky Blue","Dark Blue"))

pdf(file=paste0(dir,"power_nomisclass_2eqtl_alpha2.5e-6.pdf"),
    width=9, height=6)
print(gg)
dev.off()

# no misclass 
# 2 eqtls
# SLIGHTLY MODIFY ORDER OF BARS

tmp <- pdat[pdat$mcpt %in% c("0%","NA"),]
tmp <- tmp[tmp$testing_n==10000,]
tmp <- tmp[tmp$causal_num==2,]
tmp$rho_label <- paste0("rho == ", tmp$rho)
tmp$he2 <- paste0("afr",tmp$he2_afr,"eur",tmp$he2_eur)
tmp$he2 <- factor(tmp$he2,levels=c('afr0.1eur0.2', 'afr0.2eur0.1', 'afr0.2eur0.2'),
                  labels=c('He2 0.1/0.2',
                           'He2 0.2/0.1',
                           'He2 0.2/0.2'))
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10),labels=c("2 eQTLs","10 eQTLs"))
tmp$op <- factor(tmp$op,levels=c(0.5,1),labels=c("OP 0.5","OP 1"))
tmp <- tmp[tmp$method %in% c("casPS, aspPS, ref PS",
                             "aspPS AFR",
                             "aspPS EUR",
                             "casPS",
                             "ADMIX_80_10g PS",
                             "ADMIX_50_5g PS",
                             "AFR PS",
                             "EUR PS"),]
tmp$plot_method <- "Level 1 (LA-Unaware)"
tmp$plot_method[tmp$method %in% c("casPS, aspPS",
                                  "casPS, aspPS, ref PS")] <- "Level 2 (CADET)"
tmp$plot_method[tmp$method %in% c("aspPS AFR",
                                  "aspPS EUR",
                                  "casPS")] <- "Level 1 (LA-Aware)"
tmp$plot_method <- factor(tmp$plot_method,levels=c("Level 1 (LA-Unaware)",
                                                   "Level 1 (LA-Aware)",
                                                   "Level 2 (CADET)"))
tmp$method <- factor(tmp$method,levels=c("casPS, aspPS, ref PS",
                                         "casPS",
                                         "aspPS AFR",
                                         "aspPS EUR",
                                         
                                         "ADMIX_80_10g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"),
                     labels=c("CADET",
                              "casPS",
                              "aspPS AFR",
                              "aspPS EUR",
                              "ADMIX_80_10g PS",
                              "ADMIX_50_5g PS",
                              "AFR PS",
                              "EUR PS"))
tmp$rho_label <- paste0("rho == ", tmp$rho)
gg <- ggplot(data=tmp, aes(x=method, y=power,fill=plot_method)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
  labs(y=paste0("Power at alpha ",alpha),title="2 eQTLs, No LA Misclassification",subtitle="Testing N 10000",fill="Method",x="")+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black"),
        legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        legend.position="bottom")+
  coord_flip()+
  facet_grid(he2~op+rho_label, labeller = labeller(rho_label = label_parsed))+
  scale_fill_manual(values=c("Orchid","Deep Sky Blue","Dark Blue"))

pdf(file=paste0(dir,"power_method_order_change_nomisclass_2eqtl_alpha2.5e-6.pdf"),
    width=9, height=6)
print(gg)
dev.off()

# 10 eqtls
tmp <- pdat[pdat$mcpt %in% c("0%","NA"),]
tmp <- tmp[tmp$testing_n==10000,]
tmp <- tmp[tmp$causal_num==10,]
tmp$he2 <- paste0("afr",tmp$he2_afr,"eur",tmp$he2_eur)
tmp$he2 <- factor(tmp$he2,levels=c('afr0.1eur0.2', 'afr0.2eur0.1', 'afr0.2eur0.2'),
                  labels=c('He2 0.1/0.2',
                           'He2 0.2/0.1',
                           'He2 0.2/0.2'))
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10),labels=c("2 eQTLs","10 eQTLs"))
tmp$op <- factor(tmp$op,levels=c(0.5,1),labels=c("OP 0.5","OP 1"))
tmp <- tmp[tmp$method %in% c("casPS, aspPS, ref PS",
                             "aspPS AFR",
                             "aspPS EUR",
                             "casPS",
                             "ADMIX_80_10g PS",
                             "ADMIX_50_5g PS",
                             "AFR PS",
                             "EUR PS"),]
tmp$plot_method <- "Level 1 (LA-Unaware)"
tmp$plot_method[tmp$method %in% c("casPS, aspPS",
                                  "casPS, aspPS, ref PS")] <- "Level 2 (CADET)"
tmp$plot_method[tmp$method %in% c("aspPS AFR",
                                  "aspPS EUR",
                                  "casPS")] <- "Level 1 (LA-Aware)"
tmp$plot_method <- factor(tmp$plot_method,levels=c("Level 1 (LA-Unaware)",
                                                   "Level 1 (LA-Aware)",
                                                   "Level 2 (CADET)"))
tmp$method <- factor(tmp$method,levels=c("casPS, aspPS, ref PS",
                                         "aspPS AFR",
                                         "aspPS EUR",
                                         "casPS",
                                         "ADMIX_80_10g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"),
                     labels=c("CADET",
                              "aspPS AFR",
                              "aspPS EUR",
                              "casPS",
                              "ADMIX_80_10g PS",
                              "ADMIX_50_5g PS",
                              "AFR PS",
                              "EUR PS"))
tmp$rho_label <- paste0("rho == ", tmp$rho)
gg <- ggplot(data=tmp, aes(x=method, y=power,fill=plot_method)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
  labs(y=paste0("Power at alpha ",alpha),title="10 eQTLs, No LA Misclassification",subtitle="Testing N 10000",fill="Method",x="")+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text = element_text(size = 9),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black"),
        legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        legend.position="bottom")+
  coord_flip()+
  facet_grid(he2~op+rho_label, labeller = labeller(rho_label = label_parsed))+
  scale_fill_manual(values=c("Orchid","Deep Sky Blue","Dark Blue"))

pdf(file=paste0(dir,"power_nomisclass_10eqtl_alpha2.5e-6.pdf"),
    width=9, height=6)
print(gg)
dev.off()

# 100 eqtls
tmp <- pdat[pdat$mcpt %in% c("0%","NA"),]
tmp <- tmp[tmp$testing_n==10000,]
tmp <- tmp[tmp$causal_num==100,]
tmp$he2 <- paste0("afr",tmp$he2_afr,"eur",tmp$he2_eur)
tmp$he2 <- factor(tmp$he2,levels=c('afr0.1eur0.2', 'afr0.2eur0.1', 'afr0.2eur0.2'),
                  labels=c('He2 0.1/0.2',
                           'He2 0.2/0.1',
                           'He2 0.2/0.2'))
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10,100),labels=c("2 eQTLs","10 eQTLs","100 eQTLs"))
tmp$op <- factor(tmp$op,levels=c(0.5,1),labels=c("OP 0.5","OP 1"))
tmp <- tmp[tmp$method %in% c("casPS, aspPS, ref PS",
                             "aspPS AFR",
                             "aspPS EUR",
                             "casPS",
                             "ADMIX_80_10g PS",
                             "ADMIX_50_5g PS",
                             "AFR PS",
                             "EUR PS"),]
tmp$plot_method <- "Level 1 ACAT: standard PS"
tmp$plot_method[tmp$method %in% c("casPS, aspPS",
                                  "casPS, aspPS, ref PS")] <- "Level 2 ACAT"
tmp$plot_method[tmp$method %in% c("aspPS AFR",
                                  "aspPS EUR",
                                  "casPS")] <- "Level 1 ACAT: LA-aware PS"
tmp$plot_method <- factor(tmp$plot_method,levels=c("Level 1 ACAT: standard PS",
                                                   "Level 1 ACAT: LA-aware PS",
                                                   "Level 2 ACAT"))
tmp$method <- factor(tmp$method,levels=c("casPS, aspPS, ref PS",
                                         "aspPS AFR",
                                         "aspPS EUR",
                                         "casPS",
                                         "ADMIX_80_10g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"),
                     labels=c("CADET",
                              "aspPS AFR",
                              "aspPS EUR",
                              "casPS",
                              "ADMIX_80_10g PS",
                              "ADMIX_50_5g PS",
                              "AFR PS",
                              "EUR PS"))
tmp$rho_label <- paste0("rho == ", tmp$rho)
gg <- ggplot(data=tmp, aes(x=method, y=power,fill=plot_method)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
  labs(y=paste0("Power at alpha ",alpha),title="100 eQTLs, No LA Misclassification",subtitle="Testing N 10000",fill="Method",x="")+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        legend.position="bottom")+
  coord_flip()+
  facet_grid(he2~op+rho_label, labeller = labeller(rho_label = label_parsed))+
  scale_fill_manual(values=c("Orchid","Deep Sky Blue","Dark Blue"))

pdf(file=paste0(dir,"power_nomisclass_100eqtl_alpha2.5e-6.pdf"),
    width=9, height=6)
print(gg)
dev.off()

# 2 eqtls 10% misclassification

tmp <- pdat[pdat$mcpt %in% c("10%","NA"),]
tmp <- tmp[tmp$testing_n==10000,]
tmp <- tmp[tmp$causal_num==2,]
tmp$he2 <- paste0("afr",tmp$he2_afr,"eur",tmp$he2_eur)
tmp$he2 <- factor(tmp$he2,levels=c('afr0.1eur0.2', 'afr0.2eur0.1', 'afr0.2eur0.2'),
                  labels=c('He2 0.1/0.2',
                           'He2 0.2/0.1',
                           'He2 0.2/0.2'))
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10),labels=c("2 eQTLs","10 eQTLs"))
tmp$op <- factor(tmp$op,levels=c(0.5,1),labels=c("OP 0.5","OP 1"))
tmp <- tmp[tmp$method %in% c("casPS, aspPS, ref PS",
                             "aspPS AFR",
                             "aspPS EUR",
                             "casPS",
                             "ADMIX_80_10g PS",
                             "ADMIX_50_5g PS",
                             "AFR PS",
                             "EUR PS"),]
tmp$plot_method <- "Level 1 ACAT: standard PS"
tmp$plot_method[tmp$method %in% c("casPS, aspPS",
                                  "casPS, aspPS, ref PS")] <- "Level 2 ACAT"
tmp$plot_method[tmp$method %in% c("aspPS AFR",
                                  "aspPS EUR",
                                  "casPS")] <- "Level 1 ACAT: LA-aware PS"
tmp$plot_method <- factor(tmp$plot_method,levels=c("Level 1 ACAT: standard PS",
                                                   "Level 1 ACAT: LA-aware PS",
                                                   "Level 2 ACAT"))
tmp$method <- factor(tmp$method,levels=c("casPS, aspPS, ref PS",
                                         "aspPS AFR",
                                         "aspPS EUR",
                                         "casPS",
                                         "ADMIX_80_10g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"),
                     labels=c("CADET",
                              "aspPS AFR",
                              "aspPS EUR",
                              "casPS",
                              "ADMIX_80_10g PS",
                              "ADMIX_50_5g PS",
                              "AFR PS",
                              "EUR PS"))
tmp$rho_label <- paste0("rho == ", tmp$rho)
gg <- ggplot(data=tmp, aes(x=method, y=power,fill=plot_method)) +
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+
  labs(y=paste0("Power at alpha ",alpha),title="2 eQTLs, 10% LA Misclassification",subtitle="Testing N 10000",fill="Method",x="")+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        legend.position="bottom")+
  coord_flip()+
  facet_grid(he2~op+rho_label, labeller = labeller(rho_label = label_parsed))+
  scale_fill_manual(values=c("Orchid","Deep Sky Blue","Dark Blue"))

pdf(file=paste0(dir,"power_10misclass_2eqtl_alpha2.5e-6.pdf"),
    width=9, height=6)
print(gg)
dev.off()

# 10 eqtls 10% misclass

tmp <- pdat[pdat$mcpt %in% c("10%","NA"),]
tmp <- tmp[tmp$testing_n==10000,]
tmp <- tmp[tmp$causal_num==10,]
tmp$he2 <- paste0("afr",tmp$he2_afr,"eur",tmp$he2_eur)
tmp$he2 <- factor(tmp$he2,levels=c('afr0.1eur0.2', 'afr0.2eur0.1', 'afr0.2eur0.2'),
                  labels=c('He2 0.1/0.2',
                           'He2 0.2/0.1',
                           'He2 0.2/0.2'))
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10),labels=c("2 eQTLs","10 eQTLs"))
tmp$op <- factor(tmp$op,levels=c(0.5,1),labels=c("OP 0.5","OP 1"))
tmp <- tmp[tmp$method %in% c("casPS, aspPS, ref PS",
                             "aspPS AFR",
                             "aspPS EUR",
                             "casPS",
                             "ADMIX_80_10g PS",
                             "ADMIX_50_5g PS",
                             "AFR PS",
                             "EUR PS"),]
tmp$plot_method <- "Level 1 ACAT: standard PS"
tmp$plot_method[tmp$method %in% c("casPS, aspPS",
                                  "casPS, aspPS, ref PS")] <- "Level 2 ACAT"
tmp$plot_method[tmp$method %in% c("aspPS AFR",
                                  "aspPS EUR",
                                  "casPS")] <- "Level 1 ACAT: LA-aware PS"
tmp$plot_method <- factor(tmp$plot_method,levels=c("Level 1 ACAT: standard PS",
                                                   "Level 1 ACAT: LA-aware PS",
                                                   "Level 2 ACAT"))
tmp$method <- factor(tmp$method,levels=c("casPS, aspPS, ref PS",
                                         "aspPS AFR",
                                         "aspPS EUR",
                                         "casPS",
                                         "ADMIX_80_10g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"),
                     labels=c("CADET",
                              "aspPS AFR",
                              "aspPS EUR",
                              "casPS",
                              "ADMIX_80_10g PS",
                              "ADMIX_50_5g PS",
                              "AFR PS",
                              "EUR PS"))
tmp$rho_label <- paste0("rho == ", tmp$rho)
gg <- ggplot(data=tmp, aes(x=method, y=power,fill=plot_method)) +
  geom_bar(stat="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y=paste0("Power at alpha ",alpha),title="10 eQTLs, 10% LA Misclassification",subtitle="Testing N 10000",fill="Method",x="")+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        legend.position="bottom")+
  coord_flip()+
  facet_grid(he2~op+rho_label, labeller = labeller(rho_label = label_parsed))+
  scale_fill_manual(values=c("Orchid","Deep Sky Blue","Dark Blue"))

pdf(file=paste0(dir,"power_10misclass_10eqtl_alpha2.5e-6.pdf"),
    width=9, height=6)
print(gg)
dev.off()

# 100 eQTLs 10% misclass

tmp <- pdat[pdat$mcpt %in% c("10%","NA"),]
tmp <- tmp[tmp$testing_n==10000,]
tmp <- tmp[tmp$causal_num==100,]
tmp$he2 <- paste0("afr",tmp$he2_afr,"eur",tmp$he2_eur)
tmp$he2 <- factor(tmp$he2,levels=c('afr0.1eur0.2', 'afr0.2eur0.1', 'afr0.2eur0.2'),
                  labels=c('He2 0.1/0.2',
                           'He2 0.2/0.1',
                           'He2 0.2/0.2'))
tmp$causal_num <- factor(tmp$causal_num,levels=c(2,10,100),labels=c("2 eQTLs","10 eQTLs","100 eQTLs"))
tmp$op <- factor(tmp$op,levels=c(0.5,1),labels=c("OP 0.5","OP 1"))
tmp <- tmp[tmp$method %in% c("casPS, aspPS, ref PS",
                             "aspPS AFR",
                             "aspPS EUR",
                             "casPS",
                             "ADMIX_80_10g PS",
                             "ADMIX_50_5g PS",
                             "AFR PS",
                             "EUR PS"),]
tmp$plot_method <- "Level 1 ACAT: standard PS"
tmp$plot_method[tmp$method %in% c("casPS, aspPS",
                                  "casPS, aspPS, ref PS")] <- "Level 2 ACAT"
tmp$plot_method[tmp$method %in% c("aspPS AFR",
                                  "aspPS EUR",
                                  "casPS")] <- "Level 1 ACAT: LA-aware PS"
tmp$plot_method <- factor(tmp$plot_method,levels=c("Level 1 ACAT: standard PS",
                                                   "Level 1 ACAT: LA-aware PS",
                                                   "Level 2 ACAT"))
tmp$method <- factor(tmp$method,levels=c("casPS, aspPS, ref PS",
                                         "aspPS AFR",
                                         "aspPS EUR",
                                         "casPS",
                                         "ADMIX_80_10g PS",
                                         "ADMIX_50_5g PS",
                                         "AFR PS",
                                         "EUR PS"),
                     labels=c("CADET",
                              "aspPS AFR",
                              "aspPS EUR",
                              "casPS",
                              "ADMIX_80_10g PS",
                              "ADMIX_50_5g PS",
                              "AFR PS",
                              "EUR PS"))
tmp$rho_label <- paste0("rho == ", tmp$rho)
gg <- ggplot(data=tmp, aes(x=method, y=power,fill=plot_method)) +
  geom_bar(stat="identity")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y=paste0("Power at alpha ",alpha),title="100 eQTLs, 10% LA Misclassification",subtitle="Testing N 10000",fill="Method",x="")+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text = element_text(size = 9),
        legend.text=element_text(size=11),
        legend.title=element_text(size=11),
        legend.position="bottom")+
  coord_flip()+
  facet_grid(he2~op+rho_label, labeller = labeller(rho_label = label_parsed))+
  scale_fill_manual(values=c("Orchid","Deep Sky Blue","Dark Blue"))

pdf(file=paste0(dir,"power_10misclass_100eqtl_alpha2.5e-6.pdf"),
    width=9, height=6)
print(gg)
dev.off()



##############################################################################
# UKB results
##############################################################################

# dat_bmi <- {}
# 
# for(chr in 1:22){
#   tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass7/results/pval_BMI_chr",chr,"_log_pheno.txt"), 
#                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#   tmp$chr <- chr
#   dat_bmi <- rbind(dat_bmi,tmp)
#   rm(tmp)
# }
# head(dat_bmi)

dat <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass6/results/pval_chr",chr,"_log_pheno.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat <- rbind(dat,tmp)
  rm(tmp)
}
head(dat)

gene_list <- unique(dat$gene)
length(gene_list) # 15093

dat_AA <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass6/results/pval_AA_chr",chr,"_log_pheno.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat_AA <- rbind(dat_AA,tmp)
  rm(tmp)
}
head(dat_AA)
dat_comb <- rbind(dat,dat_AA)
length(unique(dat_AA$gene))

gene_list <- unique(dat$gene)
length(gene_list) # 15093

tapply(dat$p,INDEX=dat$`ACAT Approach`,FUN=function(p) sum(p<0.05/(29*15093)))

# aspPS                aspPS-LAstd                      casPS                casPS-LAstd 
# 82                         57                         53                         43 
# casPS-LAstd+aspPS-LAstd casPS-LAstd+aspPS-LAstd+PS                casPS+aspPS             casPS+aspPS+PS 
# 58                         82                                        83                         82 
# PS-AFR                     PS-EUR 
# 41                        43

tapply(dat$p,INDEX=dat$`ACAT Approach`,FUN=function(p) sum(p<0.05/(15093)))
tapply(dat$gene,INDEX=dat$`ACAT Approach`,FUN=function(g) length(unique(g)))


sig <- dat[dat$p<0.05/(length(gene_list)*29),]
sig <- sig[-grep("LAstd",sig$`ACAT Approach`),]
sig$gene_trait <- paste0(sig$gene,"/",sig$pheno)
dim(sig) # 384
length(unique(sig$pheno)) # 15
table(sig$pheno)[order(table(sig$pheno),decreasing = T)]
# SHBG           Total bilirubin          Direct bilirubin             Lipoprotein A          Apolipoprotein B 
# 83                        51                        38                        35                        33 
# Gamma glutamyltransferase                LDL direct      Alkaline phosphatase        C-reactive protein               Cholesterol 
# 33                                         23                        20                        20                        17 
# Urate                Cystatin C          Apolipoprotein A                 Phosphate                Creatinine 
# 12                         6                         5                         5                         3
length(unique(sig$gene_trait)) 
# 92

# upset plot!!

# Convert the tibble to a binary membership matrix
upset_data <- sig[!sig$`ACAT Approach`=="casPS",] %>%
  select(`ACAT Approach`, gene_trait) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = `ACAT Approach`, values_from = value, values_fill = 0)

upset_df <- as.data.frame(upset_data)

gg <- upset(upset_df, sets = colnames(upset_df)[-1], 
      order.by = "freq", 
      empty.intersections = "on", 
      mainbar.y.label = "Gene-Trait Pair Intersections",
      sets.x.label = "Significant Gene-Trait Pairs",
      set_size.show=T,nintersects=12,set_size.scale_max=88)

pdf(file="/Volumes/walden/Admix_TWAS_Project/manuscript/v10_post_review_1/figures/upset.pdf",
    width=6, height=6)
print(gg)
dev.off()

# upset plot with AA results

sig <- dat[dat$p<0.05/(length(gene_list)*29),]
sig <- sig[-grep("LAstd",sig$`ACAT Approach`),]
sig$gene_trait <- paste0(sig$gene,"/",sig$pheno)
dim(sig) # 384

tmp <- dat_AA[dat_AA$p<0.05/(length(gene_list)*29),]
tmp$gene_trait <- paste0(tmp$gene,"/",tmp$pheno)
sig_comb <- rbind(sig,tmp)
length(unique(sig_comb$pheno)) # 15
table(sig_comb$pheno)[order(table(sig_comb$pheno),decreasing = T)]

# SHBG           Total bilirubin          Direct bilirubin             Lipoprotein A 
# 95                        58                        44                        39 
# Gamma glutamyltransferase          Apolipoprotein B                LDL direct      Alkaline phosphatase 
# 37                        35                        26                        23 
# C-reactive protein               Cholesterol                     Urate                Cystatin C 
# 22                        19                        15                         8 
# Phosphate          Apolipoprotein A                Creatinine 
# 6                         5                         3 


upset_data <- sig_comb[!sig_comb$`ACAT Approach`=="casPS",] %>%
  select(`ACAT Approach`, gene_trait) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = `ACAT Approach`, values_from = value, values_fill = 0)

upset_df <- as.data.frame(upset_data)

gg <- upset(upset_df, sets = colnames(upset_df)[-1], 
            order.by = "freq", 
            empty.intersections = "on", 
            mainbar.y.label = "Gene-Trait Pair Intersections",
            sets.x.label = "Significant Gene-Trait Pairs",
            set_size.show=T,nintersects=16,set_size.scale_max=88)

pdf(file="/Volumes/walden/Admix_TWAS_Project/manuscript/v10_post_review_1/figures/upset_with_AA.pdf",
    width=6, height=6)
print(gg)
dev.off()

## venn diagram

venn.diagram(
  x = list(unique(sig$gene_trait[sig$`ACAT Approach`=="PS-AA"]),
           unique(sig$gene_trait[sig$`ACAT Approach`=="casPS+aspPS"]), 
           unique(sig$gene_trait[sig$`ACAT Approach`=="PS-AFR"]),
           unique(sig$gene_trait[sig$`ACAT Approach`=="PS-EUR"]),
         #  unique(sig$gene_trait[sig$`ACAT Approach`=="casPS+aspPS+PS"]),
           unique(sig$gene_trait[sig$`ACAT Approach`=="aspPS"])),
  category.names = c("PS-AA","casPS+asPS", "PS-AFR" , "PS-EUR","aspPS"),
  filename = '/Volumes/walden/Admix_TWAS_Project/manuscript/v10_post_review_1/figures/venn_diagram_all.png',
  output=TRUE,main.fontfamily = "serif",fontfamily="Arial")

unique <- unique(sig$gene_trait)
unique_trait <- unlist(lapply(strsplit(unique,"/"), `[[`, 2))
sort(table(unique_trait),decreasing = T)
# SHBG           Total bilirubin          Apolipoprotein B             Lipoprotein A          Direct bilirubin 
# 19                        11                         9                         9                         8 
# Gamma glutamyltransferase                LDL direct               Cholesterol        C-reactive protein      Alkaline phosphatase 
# 7                                          7                         6                         5                         4 
# Apolipoprotein A                     Urate                Creatinine                Cystatin C                 Phosphate 
# 2                                   2                         1                         1                         1

headfile <- {}
for(chr in 1:22){
  d <- fread(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass6/files_for_analysis/gene_anno_files/gene_anno_all_chr",chr,".txt"))
  headfile <- rbind(headfile,d)
}
# headfile <- read_delim("/Volumes/walden/data_GTEx/expression/headfile.txt", 
#                        delim = "\t", escape_double = FALSE, 
#                        trim_ws = TRUE)
# headfile$gene <- unlist(lapply(strsplit(headfile$TargetID,"[.]"), `[[`, 1))
dat <- merge(x=sig,y=headfile,by.x="gene",by.y="TargetID",all.x=T,all.y=F)

load("/Volumes/walden/Admix_TWAS_Project/UKB/pass2/files_for_analysis/UKB_adj_phenodat_norm.RData")
pheno_vec <- names(adj_dat_norm)[-1]
pheno_vec <- pheno_vec[-which(pheno_vec=="Oestradiol")]

pheno_dat <- data.frame(pheno=pheno_vec, group=NA)
pheno_dat$group[pheno_dat$pheno %in% c("Alkaline phosphatase",
                                       "Calcium",
                                       "Vitamin D",
                                       "Rheumatoid factor")] <- "Bone and joint"
pheno_dat$group[pheno_dat$pheno %in% c("Apolipoprotein A",
                                       "Apolipoprotein B",
                                       "C-reactive protein",
                                       "Cholesterol",
                                       "HDL cholesterol",
                                       "LDL direct",
                                       "Lipoprotein A",
                                       "Triglycerides")] <- "Cardiovascular"
pheno_dat$group[pheno_dat$pheno %in% c("Glucose",
                                       "Glycated haemoglobin (HbA1c)")] <- "Diabetes"
pheno_dat$group[pheno_dat$pheno %in% c("IGF-1",
                                       "SHBG",
                                       "Testosterone")] <- "Hormone"
pheno_dat$group[pheno_dat$pheno %in% c("Alanine aminotransferase",
                                       "Albumin",
                                       "Aspartate aminotransferase",
                                       "Direct bilirubin",
                                       "Gamma glutamyltransferase",
                                       "Total bilirubin")] <- "Liver"
pheno_dat$group[pheno_dat$pheno %in% c("Creatinine",
                                       "Cystatin C",
                                       "Phosphate",
                                       "Total protein",
                                       "Urate",
                                       "Urea")] <- "Renal"
dat <- merge(x=dat,y=pheno_dat,by="pheno")


gene_LA_aware <- unique(dat$gene_trait[dat$`ACAT Approach` %in% c("casPS+aspPS")])
length(gene_LA_aware) # 83
gene_std <- unique(dat$gene_trait[dat$`ACAT Approach` %in% c("PS-EUR","PS-AFR")])
length(gene_std) # 71
focus_gene <- gene_LA_aware[!gene_LA_aware %in% gene_std]
length(focus_gene) # 18

# compare to kitchen sink results
gene_kitchen <- unique(dat$gene_trait[dat$`ACAT Approach` %in% c("casPS+aspPS+PS")])
length(gene_kitchen) # 82
gene_kitchen[!gene_kitchen %in% gene_LA_aware]
# [1] "ENSG00000171786/C-reactive protein" "ENSG00000189114/Cholesterol"        "ENSG00000104856/LDL direct"        
# [4] "ENSG00000238917/SHBG"

# main body table 1
tmp <- dat[dat$gene_trait %in% focus_gene & dat$`ACAT Approach` %in%c("casPS+aspPS"),]
tmp$CHROM <- as.numeric(tmp$CHROM)
tmp$GeneStart <- as.numeric(tmp$GeneStart)
tmp <- tmp[order(tmp$GeneStart),]
tmp <- tmp[order(tmp$CHROM),]
tmp <- tmp[,c("pheno","group","gene","GeneName","GeneType","CHROM","GeneStart","p")]
tmp$CHROM <- as.character(tmp$CHROM)
tmp$GeneStart <- as.character(tmp$GeneStart)
write_xlsx(tmp, path = paste0(out_dir,"main_table1.xlsx"))

# alt main body table 1: with only independent regions
# function to find genes within 1MB on the same chromosome
within_1MB <- function(df) {
  df$GeneStart <- as.numeric(df$GeneStart)
  result <- data.frame()
  
  # Iterate over each chromosome
  for (chr in unique(df$CHROM)) {
    chr_df <- df[df$CHROM == chr, ]
    
    # Generate pairwise combinations of rows
    if(nrow(chr_df)>1){
      for (i in 1:(nrow(chr_df) - 1)) {
        for (j in (i + 1):nrow(chr_df)) {
          dist <- abs(chr_df$GeneStart[i] - chr_df$GeneStart[j])
          if (dist <= 1000000) {  # 1MB = 1,000,000 base pairs
            result <- rbind(result, cbind(chr_df[i, ], chr_df[j, ], Distance = dist))
          }
        }
      }
    }else{
      result <- rbind(result, cbind(chr_df, chr_df, Distance = 0))
    }
    
  }
  
  return(result)
}
tab1 <- within_1MB(tmp)

# casPS+aspPS res
tmp <- dat[dat$`ACAT Approach` %in%c("casPS+aspPS"),]
tmp$CHROM <- as.numeric(tmp$CHROM)
tmp$GeneStart <- as.numeric(tmp$GeneStart)
tmp <- tmp[order(tmp$GeneStart),]
tmp <- tmp[order(tmp$CHROM),]
tmp <- tmp[,c("pheno","gene","GeneName","GeneType","CHROM","GeneStart","p")]
tmp$CHROM <- as.character(tmp$CHROM)
tmp$GeneStart <- as.character(tmp$GeneStart)
write_xlsx(tmp, path = paste0(out_dir,"stable_casPS_aspPS.xlsx"))

# casPS+aspPS+PS res
tmp <- dat[dat$`ACAT Approach` %in%c("casPS+aspPS+PS"),]
tmp$CHROM <- as.numeric(tmp$CHROM)
tmp$GeneStart <- as.numeric(tmp$GeneStart)
tmp <- tmp[order(tmp$GeneStart),]
tmp <- tmp[order(tmp$CHROM),]
tmp <- tmp[,c("pheno","gene","GeneName","GeneType","CHROM","GeneStart","p")]
tmp$CHROM <- as.character(tmp$CHROM)
tmp$GeneStart <- as.character(tmp$GeneStart)
write_xlsx(tmp, path = paste0(out_dir,"stable_casPS_aspPS_PS.xlsx"))

# aspPS res
tmp <- dat[dat$`ACAT Approach` %in%c("aspPS"),]
tmp$CHROM <- as.numeric(tmp$CHROM)
tmp$GeneStart <- as.numeric(tmp$GeneStart)
tmp <- tmp[order(tmp$GeneStart),]
tmp <- tmp[order(tmp$CHROM),]
tmp <- tmp[,c("pheno","gene","GeneName","GeneType","CHROM","GeneStart","p")]
tmp$CHROM <- as.character(tmp$CHROM)
tmp$GeneStart <- as.character(tmp$GeneStart)
write_xlsx(tmp, path = paste0(out_dir,"stable_aspPS.xlsx"))

# PS EUR res
tmp <- dat[dat$`ACAT Approach` %in%c("PS-EUR"),]
tmp$CHROM <- as.numeric(tmp$CHROM)
tmp$GeneStart <- as.numeric(tmp$GeneStart)
tmp <- tmp[order(tmp$GeneStart),]
tmp <- tmp[order(tmp$CHROM),]
tmp <- tmp[,c("pheno","gene","GeneName","GeneType","CHROM","GeneStart","p")]
tmp$CHROM <- as.character(tmp$CHROM)
tmp$GeneStart <- as.character(tmp$GeneStart)
write_xlsx(tmp, path = paste0(out_dir,"stable_PS_EUR.xlsx"))

# PS AFR res
tmp <- dat[dat$`ACAT Approach` %in%c("PS-AFR"),]
tmp$CHROM <- as.numeric(tmp$CHROM)
tmp$GeneStart <- as.numeric(tmp$GeneStart)
tmp <- tmp[order(tmp$GeneStart),]
tmp <- tmp[order(tmp$CHROM),]
tmp <- tmp[,c("pheno","gene","GeneName","GeneType","CHROM","GeneStart","p")]
tmp$CHROM <- as.character(tmp$CHROM)
tmp$GeneStart <- as.character(tmp$GeneStart)
write_xlsx(tmp, path = paste0(out_dir,"stable_PS_AFR.xlsx"))

#####################################
### SNP COUNTS IN SIMULATIONS
pass <- "38_review1"
dir <- paste0("/Volumes/walden/Admix_TWAS_Project/pass",pass)
# set working directory for figure output
out_dir="/Volumes/walden/Admix_TWAS_Project/manuscript/v10_post_review_1/figures/"

# read in parameter space
parameter_space <- read.table(paste0(dir,"/files_for_analysis/paramspace_full.txt"),sep="\t",header=T)
parameter_space$param_row <- 1:nrow(parameter_space)

# read in result file names
file.list <- list.files(paste0(dir,"/results/snp_counts"),recursive = T)

# process imputation R2 files
out <- lapply(file.list,FUN=read_in_snp_cnt,dir=dir)
dat <- bind_rows(out, .id = "column_label")

library(ggplot2)

dat$eQTL_SS <- NA
dat$eQTL_SS[dat$SUBDIR=="train_ADMIX_50_5g_500"] <- "Admixed 50 5g"
dat$eQTL_SS[dat$SUBDIR=="train_ADMIX_80_10g_500"] <- "Admixed 80 10g"
dat$eQTL_SS[dat$SUBDIR=="ref_AFR_500"] <- "Ref AFR"
dat$eQTL_SS[dat$SUBDIR=="ref_EUR_500"] <- "Ref EUR"

dat <- merge(dat, parameter_space,by.x="pr",by.y="param_row")

gg <- ggplot(dat, aes(x = COUNT, fill = METHOD)) +
  geom_histogram(binwidth = 100, position = "identity", alpha = 0.6) +  # Adjust binwidth as needed
  facet_grid(causal_num~eQTL_SS, scales = "fixed") +  # Separate histograms for each SUBDIR
  theme_bw() +
  labs(fill="PRS Model",x = "SNP Count", y = "Frequency", title = "Distribution of GReX Model SNP Counts by eQTL Dataset") +
  scale_fill_brewer(palette = "Set2")  # Adjust color scheme if needed

pdf(file="/Volumes/walden/Admix_TWAS_Project/manuscript/v10_post_review_1/new_figures/sfig_snp_counts_v2.pdf",
    width=8, height=4)
print(gg)
dev.off()

################################################
# UKB+PAGE results
################################################

## ALL PHENOS

dat <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass8/results/pval_newphen_chr",chr,"_raw_pheno.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp2 <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass8/results/pval_obes_chr",chr,"_raw_pheno.txt"), 
                     delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp <- rbind(tmp,tmp2)
  tmp$chr <- chr
  dat <- rbind(dat,tmp)
  rm(tmp)
}
head(dat)

gene_list <- unique(dat$gene)
length(gene_list) # 15093

dat_page <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/PAGE/pass2/results/pval_chr",chr,"_raw_pheno_PAGE_pass2.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat_page <- rbind(dat_page,tmp)
  rm(tmp)
}
head(dat_page)

# rename
dat$pheno[dat$pheno=="Diabetes"] <- "DIABET"
dat$pheno[dat$pheno=="heart_attack"] <- "HATTACK"
dat$pheno[dat$pheno=="High_blood_pressure"] <- "HIBP"
dat$pheno[dat$pheno=="Stroke"] <- "STROKE"
dat$pheno[dat$pheno=="obesity"] <- "OBESITY"

gene_list <- unique(dat$gene)
gene_list_page <- unique(dat_page$gene)
length(gene_list_page) # 14531
length(gene_list) # 15093
length(intersect(gene_list,gene_list_page)) # 14531

dat$gene_trait <- paste0(dat$gene,"/",dat$pheno)
dat_page$gene_trait <- paste0(dat_page$gene,"/",dat_page$pheno)
dat_page <- dat_page[!dat_page$p==0,] # i dont think some logit models converged

# meta-analysis via ACAT
# Merge datasets on gene_trait and ACAT Approach
merged_data <- inner_join(dat, dat_page, by = c("gene_trait", "ACAT Approach"), suffix = c("_ukb", "_page"))
library(ACAT)
# Apply ACAT to combine p-values
merged_data <- merged_data %>%
  rowwise() %>%
  mutate(
    p_meta = ACAT(c(p_ukb, p_page))  # Combining two p-values using ACAT
  ) %>%
  ungroup()

# View results
head(merged_data)
min(merged_data$p_meta)
merged_data[merged_data$p_meta<(0.05/(14531)),]

# sig after bonferroni correction in UKB
thresh <- 0.05/(length(gene_list)*5)
sig_ukb <- dat[dat$p<thresh,]
# # A tibble: 2 × 5
# `ACAT Approach`           p gene            pheno   chr
# <chr>                 <dbl> <chr>           <chr> <int>
#   1 casPS           0.000000569 ENSG00000174306 BMI      20
# 2 PS-AFR          0.000000473 ENSG00000174306 BMI      20

thresh <- 0.05/(length(gene_list))
sig_ukb <- dat[dat$p<thresh,]

# sig after bonferroni correction in PAGE
thresh <- 0.05/(length(gene_list_page)*5)
sig_page <- dat_page[dat_page$p<thresh,]
# 0 results

# sig after less strict correction in PAGE
thresh <- 0.05/(length(gene_list_page))
sig_page <- dat_page[dat_page$p<thresh,]

min(dat_page$p[dat_page$pheno=="OBESITY"])
sig_page <- dat_page[dat_page$p<5e-5,]

# look in UKB first, then PAGE
for(thresh in c(0.05,0.005,0.0005,5e-5)){
  print(thresh)
  
  sig_ukb <- dat[dat$p<thresh,]
  x <- unique(sig_ukb$gene_trait)
  sig_page <- dat_page[dat_page$gene_trait %in% x,]
  nsig <- sum(sig_page$p < (0.05/length(x)))
  if(nsig >0){
    t <- sig_page[sig_page$p< (0.05/length(x)),]
    print(t)
    print("ukb original res:")
    print(sig_ukb[sig_ukb$gene_trait %in% t$gene_trait,])
  }
}

dat[dat$gene_trait=="ENSG00000156858/BMI" & dat$p <5e-05,]

# look in UKB first ONLY CADET sig, then PAGE
for(thresh in c(0.05,0.005,0.0005,5e-5,5e-6)){
  print(thresh)
  
  sig_ukb <- dat[dat$p<thresh & dat$`ACAT Approach`=="casPS+aspPS+PS" | dat$p<thresh & dat$`ACAT Approach`=="casPS-LAstd+aspPS-LAstd+PS",]
  x <- unique(sig_ukb$gene_trait)
  sig_page <- dat_page[dat_page$gene_trait %in% x,]
  nsig <- sum(sig_page$p < (0.05/length(x)))
  if(nsig >0){
    t <- sig_page[sig_page$p< (0.05/length(x)),]
    print(t)
    print("ukb original res:")
    print(sig_ukb[sig_ukb$gene_trait %in% t$gene_trait,])
  }
}

# look in PAGE first, then UKB
for(thresh in c(0.05,0.005,0.0005,5e-5,5e-6)){
  print(thresh)
  
  sig_page <- dat_page[dat_page$p<thresh,]
  x <- unique(sig_page$gene_trait)
  sig_ukb <- dat[dat$gene_trait %in% x,]
  nsig <- sum(sig_ukb$p < (0.05/length(x)))
  if(nsig >0){
    t <- sig_ukb[sig_ukb$p< (0.05/length(x)),]
    print(t)
    print("orig PAGE:")
    print(sig_page[sig_page$gene_trait %in% t$gene_trait,])
    
  }
}


# just counts by approach in PAGE

tmp2 <- dat_page[-grep("LAstd",dat_page$`ACAT Approach`),]
for(thresh in c(0.05,0.005,0.0005,5e-5,5e-6)){
  print(thresh)
  for(pheno in unique(dat_page$pheno)){
    print(pheno)
    tmp <- tmp2[tmp2$p<thresh,]
    tmp <- tmp[tmp$pheno==pheno,]
    table(tmp$`ACAT Approach`)
  }
}

# Filter out "LAstd" from the ACAT Approach
tmp2 <- dat_page[!grepl("LAstd", dat_page$`ACAT Approach`),]
# Define p-value thresholds
p_thresholds <- c(0.05, 0.005, 0.0005, 5e-5, 5e-6)

# Summarize data
summary_data <- expand.grid(
  pheno = unique(tmp2$pheno),
  threshold = p_thresholds,
  `ACAT Approach` = unique(tmp2$`ACAT Approach`)
) %>%
  rowwise() %>%
  mutate(count = sum(tmp2$p < threshold & tmp2$pheno == pheno & tmp2$`ACAT Approach` == `ACAT Approach`))
summary_data$`ACAT Approach` <- factor(summary_data$`ACAT Approach`,
                                       levels=c("PS-AFR","PS-EUR","aspPS","casPS","casPS+aspPS","casPS+aspPS+PS"))

# Plot results
ggplot(summary_data, aes(x = pheno, y = count, fill = `ACAT Approach`)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ threshold, scales = "free_y") +
  #scale_y_log10() +  # Use log scale if counts vary widely
  labs(title = "Significant Associations by Phenotype and P-Value Threshold",
       x = "Phenotype", y = "Number of Significant Associations") +
  theme_minimal()

# Filter out "LAstd" from the ACAT Approach
tmp2 <- dat_page[!grepl("LAstd", dat_page$`ACAT Approach`),]

# Define p-value thresholds
p_thresholds <- c(0.05, 0.005, 0.0005, 5e-5, 5e-6)

# Summarize data (aggregate across all phenotypes)
summary_data <- expand.grid(
  threshold = p_thresholds,
  `ACAT Approach` = unique(tmp2$`ACAT Approach`)
) %>%
  rowwise() %>%
  mutate(count = sum(tmp2$p < threshold & tmp2$`ACAT Approach` == `ACAT Approach`))

summary_data$`ACAT Approach` <- factor(summary_data$`ACAT Approach`,
                                       levels = c("PS-AFR", "PS-EUR", "aspPS", "casPS", "casPS+aspPS", "casPS+aspPS+PS"))

# Plot results
ggplot(summary_data, aes(x = threshold, y = count, fill = `ACAT Approach`)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10() +  # Log scale for p-value thresholds (optional)
  labs(title = "Total Significant Associations Across All Phenotypes",
       x = "P-Value Threshold", y = "Number of Significant Associations") +
  theme_minimal()

## JUST BMI

dat <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass8/results/pval_newphen_chr",chr,"_raw_pheno.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat <- rbind(dat,tmp)
  rm(tmp)
}
head(dat)

gene_list <- unique(dat$gene)
length(gene_list) # 15093

dat_page <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/PAGE/pass1/results/pval_chr",chr,"_raw_pheno_PAGE.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat_page <- rbind(dat_page,tmp)
  rm(tmp)
}
head(dat_page)

dat <- dat[dat$pheno=="BMI",]
dat_page <- dat_page[dat_page$pheno=="BMI",]

gene_list <- unique(dat$gene)
gene_list_page <- unique(dat_page$gene)
length(gene_list_page) # 14531
length(gene_list) # 15093
length(intersect(gene_list,gene_list_page)) # 14531

dat$gene_trait <- paste0(dat$gene,"/",dat$pheno)
dat_page$gene_trait <- paste0(dat_page$gene,"/",dat_page$pheno)
dat_page <- dat_page[!dat_page$p==0,] # i dont think some logit models converged

# sig after bonferroni correction in UKB
thresh <- 0.05/(length(gene_list)*1)
sig_ukb <- dat[dat$p<thresh,]
# `ACAT Approach`                      p gene            pheno   chr gene_trait         
# <chr>                            <dbl> <chr>           <chr> <int> <chr>              
#   1 casPS                      0.000000569 ENSG00000174306 BMI      20 ENSG00000174306/BMI
# 2 aspPS                      0.00000258  ENSG00000174306 BMI      20 ENSG00000174306/BMI
# 3 PS-AFR                     0.000000473 ENSG00000174306 BMI      20 ENSG00000174306/BMI
# 4 casPS+aspPS                0.00000119  ENSG00000174306 BMI      20 ENSG00000174306/BMI
# 5 casPS-LAstd+aspPS-LAstd+PS 0.00000223  ENSG00000174306 BMI      20 ENSG00000174306/BMI
# 6 casPS+aspPS+PS             0.00000108  ENSG00000174306 BMI      20 ENSG00000174306/BMI
# 7 PS-AFR                     0.00000142  ENSG00000132793 BMI      20 ENSG00000132793/BMI

# sig after bonferroni correction in PAGE
thresh <- 0.05/(length(gene_list_page)*1)
sig_page <- dat_page[dat_page$p<thresh,]
# 0 results

# look in UKB first, then PAGE
for(thresh in c(0.05,0.005,0.0005,5e-5)){
  print(thresh)
  
  sig_ukb <- dat[dat$p<thresh,]
  x <- unique(sig_ukb$gene_trait)
  sig_page <- dat_page[dat_page$gene_trait %in% x,]
  nsig <- sum(sig_page$p < (0.05/length(x)))
  if(nsig >0){
    print(sig_page[sig_page$p< (0.05/length(x)),])
  }
}

dat[dat$gene_trait=="ENSG00000156858/BMI" & dat$p <5e-05,]

# look in UKB first ONLY CADET sig, then PAGE
for(thresh in c(0.05,0.005,0.0005,5e-5,5e-6)){
  print(thresh)
  
  sig_ukb <- dat[dat$p<thresh & dat$`ACAT Approach`=="casPS+aspPS+PS" | dat$p<thresh & dat$`ACAT Approach`=="casPS-LAstd+aspPS-LAstd+PS",]
  x <- unique(sig_ukb$gene_trait)
  sig_page <- dat_page[dat_page$gene_trait %in% x,]
  nsig <- sum(sig_page$p < (0.05/length(x)))
  if(nsig >0){
    print(sig_page[sig_page$p< (0.05/length(x)),])
  }
}

# look in PAGE first, then UKB
for(thresh in c(0.05,0.005,0.0005,5e-5,5e-6)){
  print(thresh)
  
  sig_page <- dat_page[dat_page$p<thresh,]
  x <- unique(sig_page$gene_trait)
  sig_ukb <- dat[dat$gene_trait %in% x,]
  nsig <- sum(sig_ukb$p < (0.05/length(x)))
  if(nsig >0){
    print(sig_ukb[sig_ukb$p< (0.05/length(x)),])
  }
}

# just counts by approach in PAGE

tmp2 <- dat_page[-grep("LAstd",dat_page$`ACAT Approach`),]
for(thresh in c(0.05,0.005,0.0005,5e-5,5e-6)){
  print(thresh)
  for(pheno in unique(dat_page$pheno)){
    print(pheno)
    tmp <- tmp2[tmp2$p<thresh,]
    tmp <- tmp[tmp$pheno==pheno,]
    print(table(tmp$`ACAT Approach`)[order(table(tmp$`ACAT Approach`),decreasing = T)])
  }
}


###########################################
# corr between approaches in UKB
###########################################

files <- list.files("/Volumes/walden/Admix_TWAS_Project/UKB/pass6/results/grex_corr/")
res <- {}
for(i in 1:length(files)){
   load(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass6/results/grex_corr/",files[i]))
  res <- rbind(res,res_formal)
}
head(res)  
res$Variable1Num <- as.numeric(res$Variable1)
res$Variable2Num <- as.numeric(res$Variable2)

summ <- res %>%
  mutate(pair = paste(pmin(Variable1Num, Variable2Num), pmax(Variable1Num, Variable2Num), sep = "-")) %>%
  group_by(pair, model) %>%
  summarise(
    mean_correlation = mean(Correlation, na.rm = TRUE),
    sd_correlation = sd(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

summ <- summ[summ$mean_correlation<1,]
summ$Variable1 <- as.numeric(substr(summ$pair,1,1))
summ$Variable2 <- as.numeric(substr(summ$pair,3,3))
summ$Variable1Chr <- levels(res$Variable1)[summ$Variable1]
summ$Variable2Chr <- levels(res$Variable2)[summ$Variable2]

summ$Variable1Chr <- factor(summ$Variable1Chr,
                            levels=c("grex_comb_lai_std",
                                     "grex_afr_anc_std",
                                     "grex_eur_anc_std",
                                     "grex_comb_lai",
                                     "grex_afr",
                                     "grex_eur",
                                     "grex_std_afr_ss",
                                     "grex_std_eur_ss"),
                            labels=c("casPS-LAstd",
                                     "asPS-AFR-LAstd",
                                     "asPS-EUR-LAstd",
                                     "casPS",
                                     "asPS-AFR",
                                     "asPS-EUR",
                                     "PS-AFR",
                                     "PS-EUR"))
summ$Variable2Chr <- factor(summ$Variable2Chr,
                            levels=c("grex_comb_lai_std",
                                     "grex_afr_anc_std",
                                     "grex_eur_anc_std",
                                     "grex_comb_lai",
                                     "grex_afr",
                                     "grex_eur",
                                     "grex_std_afr_ss",
                                     "grex_std_eur_ss"),
                            labels=c("casPS-LAstd",
                                     "asPS-AFR-LAstd",
                                     "asPS-EUR-LAstd",
                                     "casPS",
                                     "asPS-AFR",
                                     "asPS-EUR",
                                     "PS-AFR",
                                     "PS-EUR"))

summ_wide <- summ %>%
  pivot_wider(
    names_from = model, 
    values_from = c(mean_correlation, sd_correlation)
  )

# Print the formatted table
summ_wide <- summ_wide[-grep("LAstd",summ_wide$Variable1Chr),]
summ_wide <- summ_wide[-grep("LAstd",summ_wide$Variable2Chr),]
names(summ_wide)[1:5] <- c("pair","var1","var2","Approach1","Approach2")
to_print <- summ_wide[,c(4,5,6,9,7,10,8,11)]

write_xlsx(to_print, path = paste0("/Volumes/walden/Admix_TWAS_Project/manuscript/v10_post_review_1/supp_tables/grex_cor_UKB.xlsx"))

#################################
# PAGE new phenos with  lepik
#################################
#ukb 
dat_lepik <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass10/results/pval_chr",chr,"_new_pheno_lepik.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat_lepik <- rbind(dat_lepik,tmp)
  rm(tmp)
}
table(dat_lepik$pheno)
head(dat_lepik)
dat_lepik$gene_trait <- paste0(dat_lepik$gene,"/",dat_lepik$pheno)
length(unique(dat_lepik$gene))
sum(dat_lepik$p<(0.05/(6*15093)))
dat_lepik[dat_lepik$p<(0.05/15093),]
look_for <- unique(dat_lepik$gene_trait[dat_lepik$p<(0.05/15093)])
dat_lepik[dat_lepik$p<(0.05/(6*15093)),]
dat_lepik[dat_lepik$gene_trait=="ENSG00000250696/heart_attack",]

dat_lepik <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/PAGE/pass3/results/pval_chr",chr,"_raw_pheno_PAGE_pass3_lepik.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat_lepik <- rbind(dat_lepik,tmp)
  rm(tmp)
}
head(dat_lepik)
length(unique(dat_lepik$gene))
sum(dat_lepik$p<(0.05/(11093)))

dat <- dat_lepik[order(dat_lepik$p,decreasing=F),]
head(dat,10)
dat$gene_trait <- paste0(dat$gene,"/",dat$pheno)
dat[dat$gene_trait%in%look_for,]

#################################
# UKB with GTEX vs UKB with lepik
#################################

dat <- {}
for(chr in 1:22){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass6/results/pval_chr",chr,"_log_pheno.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat <- rbind(dat,tmp)
  rm(tmp)
}
head(dat)

dat_lepik <- {}
for(chr in 1:21){
  tmp <- read_delim(paste0("/Volumes/walden/Admix_TWAS_Project/UKB/pass9/results/pval_chr",chr,"_log_pheno_lepik.txt"), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  tmp$chr <- chr
  dat_lepik <- rbind(dat_lepik,tmp)
  rm(tmp)
}
head(dat)

gene_list <- unique(dat$gene)
gene_list_lepik <- unique(dat_lepik$gene)

sig <- dat[dat$p<(0.05/(length(gene_list)*29)),]

#dat <- dat[dat$gene %in% gene_list_lepik,]
dat_lepik <- dat_lepik[dat_lepik$gene %in% gene_list,]
length(unique(dat$gene)) # 15093
length(unique(dat_lepik$gene)) # 12388

sig$gene_trait <- paste0(sig$gene,"/",sig$pheno)
sig <- sig[-grep("LAstd",sig$`ACAT Approach`),]

sig_lepik <- dat_lepik[dat_lepik$p<(0.05/(length(unique(dat$gene))*29)),]
sig_lepik$gene_trait <- paste0(sig_lepik$gene,"/",sig_lepik$pheno)
sig_lepik <- sig_lepik[-grep("LAstd",sig_lepik$`ACAT Approach`),]

unique_genetraits <- unique(sig$gene_trait)
unique_genetraits_lepik <- unique(sig_lepik$gene_trait)
length(intersect(unique_genetraits,
                 unique_genetraits_lepik))

unique_genetraits <- unique(sig$gene_trait[sig$`ACAT Approach`=="casPS+aspPS+PS"])
length(unique_genetraits) # 54
unique_genetraits_lepik <- unique(sig_lepik$gene_trait[sig_lepik$`ACAT Approach`=="casPS+aspPS+PS"])
length(unique_genetraits_lepik) # 63
length(intersect(unique_genetraits,
                 unique_genetraits_lepik)) # 35

olap <- sig
olap$gene_in_lepik <- F
olap$gene_in_lepik[olap$gene %in% gene_list_lepik] <- T
olap$lepik_gt_match <- F
olap$lepik_gt_match[olap$gene_trait %in% sig_lepik$gene_trait] <- T
olap$lepik_gt_model_match <- F
olap$gt_model <- paste0(olap$gene_trait,"/",olap$`ACAT Approach`)
sig_lepik$gt_model <- paste0(sig_lepik$gene_trait,"/",sig_lepik$`ACAT Approach`)
olap$lepik_gt_model_match[olap$gt_model %in% sig_lepik$gt_model] <- T

# number of original sig g-t associations by CADET for genes that are in lepik
length(unique(olap$gene_trait[olap$gene_in_lepik==T & olap$`ACAT Approach`=="casPS+aspPS+PS"]))
# 54

# number of these that were also identified as significant in lepik by CADET
length(unique(olap$gene_trait[olap$gene_in_lepik==T & olap$`ACAT Approach`=="casPS+aspPS+PS" & 
                                olap$lepik_gt_model_match==T]))
#35


olap <- olap[,-c(6,10)]

write_xlsx(olap, path = paste0("/Volumes/walden/Admix_TWAS_Project/manuscript/v10_post_review_1/supp_tables/lepik_res.xlsx"))

