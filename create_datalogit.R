library(mstate)
library(ggplot2)
library(cowplot)

# Load required functions
files <- list.files("../function", full.names = TRUE)
invisible(lapply(files, source))

# Load ebmt data
data("ebmt4")
ebmt=ebmt4
# Load results from the core run
result_elf<-readRDS("ebmt results/result_elf.rda")
# Extended long data from last iteration
cox_elf=read.table("ebmt results/cox_elf.csv",sep=",",head=TRUE)
logit_elf=read.table("ebmt results/logit_elf.csv",sep=",",head=TRUE)


#### Histogram for Prob of being cured from logistic analysis
data_cure=cox_elf[,c("id","pi")]
data_cure=unique(data_cure)
hist(data_cure$pi,breaks=20,xlim=c(0,1),ylim=c(0,300),main="",xlab="Probability of being cured after transplant")

#####prediction for pt id 1238
pt_data=cox_elf[which(cox_elf$id=="1238"),]

tmat0=transMat(x = list(c(2, 3, 5, 6), c(4, 5, 6), c(4, 5, 6), c(5, 6),
                        c(6), c()), names = c("Transplant", "Recovery","Adverse Event","Adverse Event\nand Recovery",
                                              "Relapse","Death"))

data_long0= msprep(data = ebmt, trans = tmat0, time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                   status = c(NA, "rec.s", "ae.s", "recae.s","rel.s", "srv.s"), keep = c("match", "proph", "year", "agecl"))

# List of time
time_list=unique(c(data_long0$Tstart,data_long0$Tstop))
time_list=sort(time_list)

# Make msfit object when cure=0 noncured
data_cov=pt_data[1,]
data_cov$cure12=0

tran_list=unique(sort(cox_elf$trans))
tran_c=tran_list[tran_list!=as.integer(tran_list)]

curecol=c("cure12","cure13","cure16","cure24","cure26","cure34","cure36","cure46")
haz_nc=make_hazlist(cox_full=result_elf$fit_cox,tran.index=1,data_cov=data_cov,time_list=time_list,tran_c=tran_c,curecol=curecol)

for (i in c(2:13)){
  haz_nc=rbind(haz_nc,make_hazlist(cox_full=result_elf$fit_cox,tran.index=i,data_cov=data_cov,time_list = time_list,tran_c=tran_c,curecol=curecol))
}

msf_nc=list(Haz=haz_nc,trans=tmat0)
attr(msf_nc,"class")="msfit"
pre_nc=probtrans(msf_nc,predt=0,variance=FALSE)

#######Make msfit object when cure=1 cured########
tran.0=tran_list[!tran_list %in% c(tran_c,round(tran_c))]
tran_c_full=sort(c(tran.0,tran_c))
haz_c=make_hazlist(cox_full=result_elf$fit_cox,tran.index=1.2,data_cov=data_cov,time_list=time_list,tran_c=tran_c,curecol=curecol)

for (i in tran_c_full[-1]){
  if(i %in% tran.0){
    list0=data.frame(time=time_list,Haz=0,trans=i)
    haz_c=rbind(haz_c,list0)
  } else {
    haz_c=rbind(haz_c,make_hazlist(cox_full=result_elf$fit_cox,tran.index=i,data_cov=data_cov,time_list = time_list,tran_c=tran_c,curecol=curecol))
  }
}

msf_c=list(Haz=haz_c,trans=tmat0)
attr(msf_c,"class")="msfit"
pre_c=probtrans(msf_c,predt=0,variance=FALSE)


############################
#when only baseline is known
############################
P_c=predict(result_elf$fit_logit, newdata=data_cov,
             type="response")
pre_nc=probtrans(msf_nc,predt=0,variance=FALSE)

pre_c=probtrans(msf_c,predt=0,variance=FALSE)

pre_com_day0=makepre_com(pre_nc=pre_nc,pre_c=pre_c,P_c=P_c)

ord=c(5,6,4,3,2,1)
plot(pre_com_day0,ord=ord,use.ggplot = TRUE,xlim=c(0,3650))+
  labs(x = "Days since transplant", y = "Probability" ,caption="(a) Prediction from State 1 (transplant) at day 0 \n"
       ) +
  scale_fill_manual(
    values = c("Transplant"="#5e4fa2", "Recovery"="#3288bd","Adverse Event"="#66c2a5",
               "Adverse Event\nand Recovery"="#fdae61",
               "Relapse"="#d53e4f","Death"="#fee08b")
  )+theme(  
    axis.title = element_text(size = 16),        
    axis.text = element_text(size = 14) ,       
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(color = c("red", rep("black", 3))),
    axis.text.y = element_text(color = "black"),
    plot.caption = element_text(hjust = 0.5, size = 16)
  )

################################################################################
#when we know this patient had reached state 3 at day 93 and remained till day 100
################################################################################
pt_data_wide=ebmt[which(ebmt$id=="1238"),]
####prepare history data till day 100
pt_data_day100=pt_data
pt_data_day100=pt_data_day100[which(pt_data_day100$from %in% c(1,3)),]
pt_data_day100[which(pt_data_day100$from==3),"Tstop"]=100
pt_data_day100[which(pt_data_day100$from==3),"time"]=100-93
pt_data_day100[which(pt_data_day100$from==3),"status"]=0

#repeat E step to calculate posterior probability of being cured
P_c_day100=E_step(data_longlong=pt_data_day100,fit_alpha=result_elf$fit_logit,fit_cox=result_elf$fit_cox,
         tran_state   = c("12","13","15","16","24","25","26","34","35","36","45","46","56"),
         tran_noncure = c("1","2","3","4","5","6","7","8","9","10","11","12","13"),
         tran_cured=c("1.2","2.2","4.2","5.2","7.2","8.2","10.2","12.2"),
         tran_correction=c("3","6","9","11","13"),id="id",t_cure=NA)$data_logit$weight[2]

pre_nc_day100=probtrans(msf_nc,predt=100,variance=FALSE)

pre_c_day100=probtrans(msf_c,predt=100,variance=FALSE)

pre_com_day100=makepre_com(pre_nc=pre_nc_day100,pre_c=pre_c_day100,P_c=P_c_day100)

plot(pre_com_day100,ord=ord, from=3, use.ggplot = TRUE, xlim=c(0, 3650)) +
  labs(x = "Days since transplant", y = "Probability",caption="(b) Prediction from State 3 (adverse event) at day 100 \n post-transplant"
       ) +
  scale_x_continuous(
    breaks = c(100,1000,2000,3000),
    labels = c(expression(bold(100)), "1000", "2000", "3000")
  ) +
  theme_minimal() +
  theme( 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(color = c("red", rep("black", 3))),
    axis.title = element_text(size = 16),        
    axis.text = element_text(size = 14) ,       
    axis.text.y = element_text(color = "black"),
    plot.caption = element_text(hjust = 0.5, size = 16)
  )+
  scale_fill_manual(
    values = c("Transplant"="#5e4fa2", "Recovery"="#3288bd","Adverse Event"="#66c2a5",
               "Adverse Event\nand Recovery"="#fdae61",
               "Relapse"="#d53e4f","Death"="#fee08b")
  ) 


#########################################################
#when we know this patient had reached state 4 at day 310 
##########################################################
pt_data_day310=pt_data[which(pt_data$from %in% c(1,3)),]

#repeat E step to calculate posterior probability of being cured
P_c_day310=E_step(data_longlong=pt_data_day310,fit_alpha=result_elf$fit_logit,fit_cox=result_elf$fit_cox,
                  tran_state   = c("12","13","15","16","24","25","26","34","35","36","45","46","56"),
                  tran_noncure = c("1","2","3","4","5","6","7","8","9","10","11","12","13"),
                  tran_cured=c("1.2","2.2","4.2","5.2","7.2","8.2","10.2","12.2"),
                  tran_correction=c("3","6","9","11","13"),id="id",t_cure=NA)$data_logit$weight[2]

pre_nc_day310=probtrans(msf_nc,predt=310,variance=FALSE)

pre_c_day310=probtrans(msf_c,predt=310,variance=FALSE)

pre_com_day310=makepre_com(pre_nc=pre_nc_day310,pre_c=pre_c_day310,P_c=P_c_day310)


plot(pre_com_day310, ord=ord,from=4, use.ggplot = TRUE, xlim=c(0, 3650)) +
  labs(x = "Days since transplant", y = "Probability",caption="(c) Prediction from State 4 (adverse event and recovery) at day 310 \n post-transplant"
       ) +
  scale_x_continuous(
    breaks = c(310,1000,2000,3000),
    labels = c(expression(bold(310)), "1000", "2000", "3000")
  ) +
  theme_minimal() +
  theme(    legend.position = "none",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(color = c("red", rep("black", 3))),
    axis.title = element_text(size = 16),         
    axis.text = element_text(size = 14) ,       
    axis.text.y = element_text(color = "black"),
    plot.caption = element_text(hjust = 0.5, size = 16)
  )+
  scale_fill_manual(
    values = c("Transplant"="#5e4fa2", "Recovery"="#3288bd","Adverse Event"="#66c2a5",
               "Adverse Event\nand Recovery"="#fdae61",
               "Relapse"="#d53e4f","Death"="#fee08b")
  ) 

