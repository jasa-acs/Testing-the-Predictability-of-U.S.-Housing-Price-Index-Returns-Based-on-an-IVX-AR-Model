 source("multi_function.R")

 #########################################
 ### The code below replicates Table 3 ###
 #########################################

 Pi=c(0.534,0.674,0.977,0.998)
 beta=c(0,0,0,0)
 phi_1=0.517
 mu=c(0,0,0,0,0)
 sigma=matrix(c(0.0412,0.0347,-0.0299,0.1088,-0.0001,
                0.0347,0.7664,0.1072,0.9606,-0.0008,
                -0.0299,0.1072,1.4576,0.4638,0.0006,
                0.1088,0.9606,0.4638,11.4784,0.0011,
                -0.0001,-0.0008,0.0006,0.0011,0.0001),ncol=5,byrow=T)
 psi_1=c(-0.116,-0.002,0.174,0.628)
 level=0.05
 rep=10000
 garch_1=0.04
 garch_2=0.95

 n=100
 apply(coverage(n,Pi,beta,phi_1,mu,sigma,psi_1,level, rep,garch_1,garch_2),2,mean)  # Sample size=100

 n=200
 apply(coverage(n,Pi,beta,phi_1,mu,sigma,psi_1,level, rep,garch_1,garch_2),2,mean)  # Sample size=200

 n=500
 apply(coverage(n,Pi,beta,phi_1,mu,sigma,psi_1,level, rep,garch_1,garch_2),2,mean)  # Sample size=500

 ## End of Table 3 replication here!





