 library(ggplot2)
 source("uni_function.R")

 ##########################################
 ### The code below replicates Figure 1 ###
 ##########################################

 # Setup of parameters
 n=200
 beta=0
 phi_1=0.7
 delta=0
 psi_1=0
 level=0.05 
 rep=10000
 garch_1=0.04
 garch_2=0.95

 # Values of beta
 index<-c(seq(0,0.25,length.out=25))

 powernum=length(index)
 plot1_fg1=rep(0,powernum*3)
 plot2_fg1=rep(0,powernum*3)
 plot3_fg1=rep(0,powernum*3)
 plot4_fg1=rep(0,powernum*3)
 plot5_fg1=rep(0,powernum*3)
 plot6_fg1=rep(0,powernum*3)
 dim(plot1_fg1)=c(powernum,3)
 dim(plot2_fg1)=c(powernum,3)
 dim(plot3_fg1)=c(powernum,3)
 dim(plot4_fg1)=c(powernum,3)
 dim(plot5_fg1)=c(powernum,3)
 dim(plot6_fg1)=c(powernum,3)

 # Calculation of rejection rates by IVX-AR(1) using AIC and BIC and IVX-KMS
 for(j in 1:powernum){
                
                 plot1_fg1[j,]=apply(coverage(n,0.2,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.2
                 plot2_fg1[j,]=apply(coverage(n,0.8,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.8
                 plot3_fg1[j,]=apply(coverage(n,0.95,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=0.95
                 plot4_fg1[j,]=apply(coverage(n,1,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)     # Pi=1
                 plot5_fg1[j,]=apply(coverage(n,1.005,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=1.005
                 plot6_fg1[j,]=apply(coverage(n,1.01,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.01
                    
 }


 ## Start of plotting from here!

 # Plot 1
 dat0<-stack(plot1_fg1)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 2
 dat0<-stack(plot2_fg1)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 3
 dat0<-stack(plot3_fg1)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 4
 dat0<-stack(plot4_fg1)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 5
 dat0<-stack(plot5_fg1)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 6
 dat0<-stack(plot6_fg1)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)


 ## End of replicating Figure 1 here!




 ##########################################
 ### The code below replicates Figure 2 ###
 ##########################################

 # Setup of parameters
 n=200
 beta=0
 phi_1=0.7
 delta=0.4
 psi_1=0.2
 level=0.05 
 rep=10000
 garch_1=0.04
 garch_2=0.95

 # Values of beta
 index<-c(seq(0,0.25,length.out=25))

 powernum=length(index)
 plot1_fg2=rep(0,powernum*3)
 plot2_fg2=rep(0,powernum*3)
 plot3_fg2=rep(0,powernum*3)
 plot4_fg2=rep(0,powernum*3)
 plot5_fg2=rep(0,powernum*3)
 plot6_fg2=rep(0,powernum*3)
 dim(plot1_fg2)=c(powernum,3)
 dim(plot2_fg2)=c(powernum,3)
 dim(plot3_fg2)=c(powernum,3)
 dim(plot4_fg2)=c(powernum,3)
 dim(plot5_fg2)=c(powernum,3)
 dim(plot6_fg2)=c(powernum,3)

 # Calculation of rejection rates by IVX-AR(1) using AIC and BIC and IVX-KMS
 for(j in 1:powernum){
                
                 plot1_fg2[j,]=apply(coverage(n,0.2,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.2
                 plot2_fg2[j,]=apply(coverage(n,0.8,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.8
                 plot3_fg2[j,]=apply(coverage(n,0.95,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=0.95
                 plot4_fg2[j,]=apply(coverage(n,1,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)     # Pi=1
                 plot5_fg2[j,]=apply(coverage(n,1.005,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=1.005
                 plot6_fg2[j,]=apply(coverage(n,1.01,beta+index[j]/1,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.01
                    
}


 ## Start of plotting from here!

 # Plot 1
 dat0<-stack(plot1_fg2)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 2
 dat0<-stack(plot2_fg2)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 3
 dat0<-stack(plot3_fg2)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 4
 dat0<-stack(plot4_fg2)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 5
 dat0<-stack(plot5_fg2)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 6
 dat0<-stack(plot6_fg2)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)
 

 ## End of replication Figure 2 here!




 ##########################################
 ### The code below replicates Figure 3 ###
 ##########################################

 # Setup of parameters
 n=200
 beta=0
 phi_1=0.55
 phi_2=0.18
 phi_3=-0.24
 phi_4=0.30
 phi_5=-0.05
 delta=0.4
 psi_1=0.2
 level=0.05
 rep=10000
 garch_1=0.04
 garch_2=0.95

 # Values of beta
 index<-c(seq(0,0.25,length.out=25))

 powernum=length(index)
 plot1_fg3=rep(0,powernum*3)
 plot2_fg3=rep(0,powernum*3)
 plot3_fg3=rep(0,powernum*3)
 plot4_fg3=rep(0,powernum*3)
 plot5_fg3=rep(0,powernum*3)
 plot6_fg3=rep(0,powernum*3)
 dim(plot1_fg3)=c(powernum,3)
 dim(plot2_fg3)=c(powernum,3)
 dim(plot3_fg3)=c(powernum,3)
 dim(plot4_fg3)=c(powernum,3)
 dim(plot5_fg3)=c(powernum,3)
 dim(plot6_fg3)=c(powernum,3)

 # Calculation of rejection rates by IVX-AR(1) using AIC and BIC and IVX-KMS
 for(j in 1:powernum){
                
                 plot1_fg3[j,]=apply(coverage(n,0.2,beta+index[j]/1,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.2
                 plot2_fg3[j,]=apply(coverage(n,0.8,beta+index[j]/1,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.8
                 plot3_fg3[j,]=apply(coverage(n,0.95,beta+index[j]/1,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=0.95
                 plot4_fg3[j,]=apply(coverage(n,1,beta+index[j]/1,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)     # Pi=1
                 plot5_fg3[j,]=apply(coverage(n,1.005,beta+index[j]/1,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=1.005
                 plot6_fg3[j,]=apply(coverage(n,1.01,beta+index[j]/1,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.01
                    
                     }


 ## Start of plotting from here!

 # Plot 1
 dat0<-stack(plot1_fg3)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 2
 dat0<-stack(plot2_fg3)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 3
 dat0<-stack(plot3_fg3)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 4
 dat0<-stack(plot4_fg3)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 5
 dat0<-stack(plot5_fg3)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)

#############################################################################################
 
 # Plot 6
 dat0<-stack(plot6_fg3)[,1]
 X<-rep(seq(0,0.25,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.16, y=0.3,label="Pi~' =1'",parse=T,size=10)
 

 ## End of replication Figure 3 here!







 
