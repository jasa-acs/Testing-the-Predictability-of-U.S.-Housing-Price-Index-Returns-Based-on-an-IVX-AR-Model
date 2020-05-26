 source("uni_function.R")

 # True beta
 beta=0

 # Hypothetical size
 level=0.05

 # Number of repetitions
 rep=10000

 # Parameters for GARCH(1,1)
 garch_1=0.04 ; garch_2=0.95



 #####################################################
 ## The parameters below fit for Panel 1 in Table 1 ##
 #####################################################

 # Copula parameter
 delta=0

 # AR(1) coefficient for e_t
 psi_1=0

 # AR(1) coefficient for u_t and sample size
 phi_1=-0.9 ; n=100  

 apply(coverage(n,1.01,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.01
 apply(coverage(n,1.005,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=1.005
 apply(coverage(n,1,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)     # Pi=1
 apply(coverage(n,0.95,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=0.95
 apply(coverage(n,0.8,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.8
 apply(coverage(n,0.2,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.2

 # To reproduce Panel 1 of Table 1, one could adjust phi_1 and n, 
 # where phi_1 is in {-0.9,-0.7,0,0.7,0.9} and n is in {100,200,500}.





 #####################################################
 ## The parameters below fit for Panel 2 in Table 1 ##
 #####################################################

 # Copula parameter
 delta=0.4

 # AR(1) coefficient for e_t
 psi_1=0.2

 # AR(1) coefficient for u_t and sample size
 phi_1=-0.9 ; n=100  

 apply(coverage(n,1.01,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.01
 apply(coverage(n,1.005,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean) # Pi=1.005
 apply(coverage(n,1,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)     # Pi=1
 apply(coverage(n,0.95,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=0.95
 apply(coverage(n,0.8,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.8
 apply(coverage(n,0.2,beta,phi_1,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.2

 # To reproduce Panel 2 of Table 1, one could adjust phi_1 and n, 
 # where phi_1 is in {-0.9,-0.7,0,0.7,0.9} and n is in {100,200,500}.


 # Replication of Table 1 ends here!







 #############################################################
 ## The parameters below fit for the five groups in Table 2 ##
 #############################################################

 # Sample size
 n=100 # Could also be 200 and 500

 # True beta
 beta=0

 # Copula parameter
 delta=0.4

 # AR(1) coefficient for e_t
 psi_1=0.2

 # Hypothetical size
 level=0.05

 # Number of repetitions
 rep=10000

 # GARCH(1,1) parameters
 garch_1=0.04 ; garch_2=0.95


 # Group i in Table 2
 phi_1=0.55 ; phi_2=0.18 ; phi_3=-0.24 ; phi_4=0.30 ; phi_5=-0.05

 apply(coverage_AR5(n,1.01,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=1.01
 apply(coverage_AR5(n,1.005,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.005
 apply(coverage_AR5(n,1,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)      # Pi=1
 apply(coverage_AR5(n,0.95,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.95
 apply(coverage_AR5(n,0.8,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.8
 apply(coverage_AR5(n,0.2,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.2

 # Group ii in Table 2
 phi_1=0.42 ; phi_2=-0.08 ; phi_3=0.33 ; phi_4=0.25 ; phi_5=-0.07

 apply(coverage_AR5(n,1.01,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=1.01
 apply(coverage_AR5(n,1.005,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.005
 apply(coverage_AR5(n,1,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)      # Pi=1
 apply(coverage_AR5(n,0.95,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.95
 apply(coverage_AR5(n,0.8,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.8
 apply(coverage_AR5(n,0.2,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.2

 # Group iii in Table 2
 phi_1=0.37 ; phi_2=-0.10 ; phi_3=0.31 ; phi_4=0.27 ; phi_5=0.04

 apply(coverage_AR5(n,1.01,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=1.01
 apply(coverage_AR5(n,1.005,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.005
 apply(coverage_AR5(n,1,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)      # Pi=1
 apply(coverage_AR5(n,0.95,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.95
 apply(coverage_AR5(n,0.8,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.8
 apply(coverage_AR5(n,0.2,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.2

 # Group iv in Table 2
 phi_1=0.24 ; phi_2=-0.03 ; phi_3=0.45 ; phi_4=0.17 ; phi_5=-0.09

 apply(coverage_AR5(n,1.01,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=1.01
 apply(coverage_AR5(n,1.005,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.005
 apply(coverage_AR5(n,1,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)      # Pi=1
 apply(coverage_AR5(n,0.95,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.95
 apply(coverage_AR5(n,0.8,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.8
 apply(coverage_AR5(n,0.2,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.2

 # Group v in Table 2
 phi_1=0.16 ; phi_2=0.08 ; phi_3=-0.25 ; phi_4=0.48 ; phi_5=-0.05

 apply(coverage_AR5(n,1.01,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=1.01
 apply(coverage_AR5(n,1.005,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)  # Pi=1.005
 apply(coverage_AR5(n,1,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)      # Pi=1
 apply(coverage_AR5(n,0.95,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)   # Pi=0.95
 apply(coverage_AR5(n,0.8,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.8
 apply(coverage_AR5(n,0.2,beta,phi_1,phi_2,phi_3,phi_4,phi_5,delta,psi_1,level, rep,garch_1,garch_2),2,mean)    # Pi=0.2

 # Replication of Table 2 ends here!

