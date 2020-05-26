 library(ggplot2)
 source("multi_function.R")

 ###########################################
 ### The code belowe replicates Figure 4 ###
 ###########################################

 n=200
 Pi=c(0.998,0.977,0.674,0.534)
 beta=c(0,0,0,0)
 phi_1=0.5165
 delta=0.4
 psi_1=c(0.6279,0.1741,-0.0016,-0.1159)
 level=0.05
 rep=10000
 garch_1=0.04
 garch_2=0.95

 index<-c(seq(0,0.25,length.out=25))

 powernum=length(index)
 plot1_fg4=rep(0,powernum*15)
 plot2_fg4=rep(0,powernum*15)
 plot3_fg4=rep(0,powernum*15)
 plot4_fg4=rep(0,powernum*15)
 dim(plot1_fg4)=c(powernum,15)
 dim(plot2_fg4)=c(powernum,15)
 dim(plot3_fg4)=c(powernum,15)
 dim(plot4_fg4)=c(powernum,15)

 # Calculation of rejection rates by IVX-AR(1) using AIC and BIC and IVX-KMS
 for(j in 1:powernum){
                
                 plot1_fg4[j,]=apply(coverage(n,Pi,c(index[j],0,0,0),phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=0.998
                 plot2_fg4[j,]=apply(coverage(n,Pi,c(0,index[j],0,0),phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=0.977
                 plot3_fg4[j,]=apply(coverage(n,Pi,c(0,0,index[j],0),phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=0.674
                 plot4_fg4[j,]=apply(coverage(n,Pi,c(0,0,0,index[j]),phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=0.534
                    
}

 ## Start of plotting from here!

 # Plot 1
 dat0<-stack(read.table("plot1_fg4.txt")[,c(5,10,15)])[,1]
 X<-rep(seq(0,0.35,length.out=25),3)
 dat1<-as.data.frame(cbind(X,dat0))
 dat1$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat1) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.3) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.25, y=0.3,label="Pi~'= 0.998'",parse=T,size=10)

#############################################################################################

 # Plot 2
 dat0<-stack(read.table("plot2_fg4.txt")[,c(5,10,15)])[,1]
 X<-rep(seq(0,0.35,length.out=25),3)
 dat95<-as.data.frame(cbind(X,dat0))
 dat95$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat95) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.25) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.25, y=0.3,label="Pi~'= 0.977'",parse=T,size=10)

#############################################################################################

 # Plot 3
 dat0<-stack(read.table("plot3_fg4.txt")[,c(5,10,15)])[,1]
 X<-rep(seq(0,0.35,length.out=25),3)
 dat8<-as.data.frame(cbind(X,dat0))
 dat8$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat8) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + theme(legend.position = "none") + 
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.25, y=0.3,label="Pi~'= 0.674'",parse=T,size=10)

#############################################################################################

 # Plot 4
 dat0<-stack(read.table("plot4_fg4.txt")[,c(5,10,15)])[,1]
 X<-rep(seq(0,0.35,length.out=25),3)
 dat2<-as.data.frame(cbind(X,dat0))
 dat2$Name<-c(rep("IVX-AR(1)-AIC",25),rep("IVX-AR(1)-BIC",25),rep("IVX-KMS",25))

 windows(7,7)
 ggplot(data=dat2) +
 theme(axis.title.x=element_blank(),
       axis.title.y=element_blank()) +
 geom_smooth(mapping=aes(x=X,y=dat0,color=Name,linetype=Name),size=1.5,se=F,span=0.4) + 
 scale_color_manual(values=c('#000000','#FF3300','blue1'))+
 scale_linetype_manual(values=c("solid","twodash","dotted")) +
 scale_size_manual(values=c(1.5, 1.5, 1.5)) + #theme(legend.position = "none") + 
 theme(legend.position = c(0.7, 0.4)) + theme(legend.text=element_text(size=20)) + 
 theme(legend.title = element_blank()) +
 geom_hline(yintercept=0.05,linetype="solid",size=0.5) + 
 annotate('text',x=0.25, y=0.3,label="Pi~'= 0.534'",parse=T,size=10)


 ## End of Figure 4 replication here!
