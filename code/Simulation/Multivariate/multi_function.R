 library(forecast)
 library(copula)
 library(doParallel)
 library(ggplot2)
 library(MASS)

########################################################################################################################
### The code below implements the IVX-KMS method based on Kostakis, A., Magdalinos, T. and Stamatogiannis, M. (2015) ###
########################################################################################################################

 IVX=function(data,xxrow,yyrow){

	xt=data[,xxrow]
 	if(class(xt)=="numeric"){
		xt=matrix(xt)
 	}
 	yt=data[,yyrow]
 	if(class(yt)=="numeric"){
		yt=matrix(yt)
 	}
 	n1=nrow(xt)
 	constant=0
 	beta=0.95
 	npow=1/3
 	cz=-1
 	ydata=yt[2:n1,]
 	if(class(ydata)=="numeric"){
		ydata=matrix(ydata)
 	}
 	xdata=xt[1:(n1-1),]
 	if(class(xdata)=="numeric"){
		xdata=matrix(xdata)
 	}
 	xdata2=xt[2:n1,]
 	if(class(xdata2)=="numeric"){
		xdata2=matrix(xdata2)
 	}
 	dxdata=xdata2-xdata

 	n=nrow(xdata)
 	ky=ncol(ydata)
 	kx=ncol(xdata)
 	p=kx+ky
 	Rnhat=matrix(data=0,nrow=kx,ncol=kx)
 	RnhatCterm=matrix(data=0,nrow=kx,ncol=1)
 	zinstr=matrix(data=0,nrow=n,ncol=kx)

 	const0=matrix(data=1,nrow=nrow(xdata),ncol=1)
 	xdataC=cbind(const0,xdata)
 	xdataC=as.matrix(xdataC)
 	Alpha1C=t(ydata)%*%xdataC%*%(solve(t(xdataC)%*%xdataC))
 	u0=ydata-xdataC%*%t(Alpha1C)
 	Alpha1=Alpha1C[,2:(kx+1)]
 	i=1
 	while(i<=kx){
		Rnhat[i,i]=(solve(t(xdata[,i])%*%(xdata[,i])))%*%t(xdata[,i])%*%(xdata2[,i])
 		i=i+1
 	}
 	xdata2=as.matrix(xdata2)
 	xdata=as.matrix(xdata)
 	ux=xdata2-xdata%*%t(Rnhat)

 	ut=cbind(u0,ux)
 	utt=t(ut)
 	sts=matrix(data=0,nrow=p,ncol=p)
 	sts3=matrix(data=0,nrow=p,ncol=p)
 	sts4=matrix(data=0,nrow=p,ncol=p)
 	sts5=matrix(data=0,nrow=p,ncol=p)
 	sts6=matrix(data=0,nrow=p,ncol=p)

 	M1=n^npow
 	M=trunc(M1)
 
	h=0
 	while(h<=M){
		t=h+1
 		while(t<=n){
			sts1=(1-(h/(M+1)))*(utt[,t])%*%t(utt[,t-h])
 			sts=sts+sts1
 			t=t+1
 		}
 		h=h+1
 	}
 	deltahat=sts/n
 
	h2=1
 	while(h2<=M){
		t2=h2+1
 		while(t2<=n){
			sts5=(1-(h2/(M+1)))*(utt[,t2])%*%t(utt[,t2-h2])
 			sts6=sts6+sts5
 			t2=t2+1
 		}
 		h2=h2+1
 	}
 	lamda=sts6/n
 	lamdax0=lamda[(ky+1):p,1:ky] 

 	t1=1
 	while(t1<=n){
		sts3=(utt[,t1])%*%t(utt[,t1])
 		sts4=sts4+sts3
 		t1=t1+1
 	}
 	sigmaut=sts4/n
 	sigmaeigenv=eigen(sigmaut)
 	if(min(sigmaeigenv$values)<0){
		print("Sigma matrix is not positive definite!")
 	}

 	zinstr1=matrix(da=0,nrow=n,ncol=kx)
 	t2=1
 	while(t2<=n){
		sts7=matrix(da=0,nrow=kx,ncol=1)
 		sts8=matrix(da=0,nrow=kx,ncol=1)
 		j=1
 		while(j<=t2){
			sts7=as.matrix(c(as.vector(((1+cz/(n^beta))^(t2-j))*t(t(dxdata[j,])))),ncol=1)
 			sts8=sts8+sts7
 			j=j+1
 		}
 		zinstr1[t2,]=t(sts8)
 		t2=t2+1
 	}
 	zinstr4=matrix(da=0,nrow=n,ncol=kx)
 	zinstr4[2:n,]=zinstr1[1:(n-1),]

 	zll=zinstr4
 	omegamatrix=deltahat+t(lamda)
 	omegaxx=omegamatrix[(ky+1):p,(ky+1):p]
 	omega00b=sigmaut[1:ky,1:ky]
 	omega0xb=sigmaut[1:ky,(ky+1):p]
 	omegafmd=omega00b-(omega0xb+t(lamdax0))%*%(solve(omegaxx))%*%t(omega0xb+t(lamdax0))

 	ydatadm=ydata-(as.matrix(apply(ydata,2,mean)%o%matrix(da=1,nrow=nrow(ydata),ncol=ky)))
 	qq=as.matrix(apply(xdata,2,mean)%o%matrix(da=1,nrow=nrow(xdata),ncol=kx))
 	DD1=matrix(NA,kx,n)
 	for(jj in 1:n){
 	       for(ii in 1:kx){
 	       	      DD1[ii,jj]=qq[(jj-1)*kx+ii,1]
 		}
 	}
 	xdatadm=xdata-t(DD1)
 	Alphatilb=(t(ydatadm)%*%zll)%*%solve(t(xdatadm)%*%zll)
 	zllbar=apply(zll,2,mean)
 	Pz=zll%*%solve(t(zll)%*%zll)%*%t(zll)
 	valphatilb=matrix(da=Alphatilb,nrow=nrow(Alphatilb)*ncol(Alphatilb),ncol=1)

 	Amatrix=matrix(da=0,nrow=ky,ncol=kx)
 	for(iy in 1:ky){
 	       for(jx in 1:kx){
 	       	      hmat2=matrix(da=0,nrow=ky,ncol=kx)
 		      hmat2[iy,jx]=1
 		      hmat1=t(matrix(da=hmat2,nrow=nrow(hmat2)*ncol(hmat2),ncol=1))
 		      Wn3indi=t(hmat1%*%valphatilb)%*%(solve(hmat1%*%(((solve(t(zll)%*%xdatadm))%x%diag(ky))%*%(t(zll)%*%zll%x%omega00b-n*(zllbar%*%t(zllbar))%x%omegafmd)%*%((solve(t(xdatadm)%*%zll))%x%diag(ky)))%*%t(hmat1)))%*%(hmat1%*%valphatilb)
 		      Amatrix[iy,jx]=Wn3indi
 		}
 	}
 	restr=ky*kx
 	Wn3=t(valphatilb)%*%(solve((((solve(t(zll)%*%xdatadm))%x%diag(ky))%*%(t(zll)%*%zll%x%omega00b-n*(zllbar%*%t(zllbar))%x%omegafmd)%*%((solve(t(xdatadm)%*%zll))%x%diag(ky)))))%*%(valphatilb)
 	Amatrixp=pchisq(Amatrix, df=1, lower.tail=FALSE)
 	Wn3p=pchisq(Wn3, df=kx, lower.tail=FALSE)
 
	coefficients=round(Alphatilb,4)
 	wald1=round(Amatrix,3)
 	pvalue=Amatrixp
 	wald=matrix(NA,1,kx)
 	for(ii in 1:kx){
 	       if(pvalue[1,ii]<=0.01){
			wald[1,ii]=paste(wald1[1,ii],"***",sep="")
 		}
 		else {if(pvalue[1,ii]<=0.05){
				wald[1,ii]=paste(wald1[1,ii],"**",sep="")
 			}
 			else {if(pvalue[1,ii]<=0.1){
					wald[1,ii]=paste(wald1[1,ii],"*",sep="")
 				}
 				else {wald[1,ii]=wald1[1,ii]
 				}
 			}
 		}
 	}

 	jointwald1=round(Wn3,3)
 	jointpvalue=Wn3p
 	if(jointpvalue[1,1]<=0.01){
 	jointwald=paste(jointwald1[1,1],"***",sep="")
 	} else{
		if(jointpvalue[1,1]<=0.05){
			jointwald=paste(jointwald1[1,1],"**",sep="")
 		} else{
			if(jointpvalue[1,1]<=0.1){
				jointwald=paste(jointwald1[1,1],"*",sep="")
 			} else{
				jointwald=jointwald1[1,1]
			} 
 		} 
 	}

 	result=matrix(NA,5,kx)
 	for(jj in 1:(kx)){
 	       result[1,jj]=coefficients[1,jj]
 	       result[2,jj]=wald[1,jj]
 	       result[3,jj]=pvalue[1,jj]
 	}
 	result[4,1]=jointwald
 	result[5,1]=jointpvalue
 	rownames(result)=c("coefficients","wald","pvalue","jointwald","jointpvalue")
 	colnames(result)=xxrow
 
	return(result)

 } 


##################################################################################################################
### The code below implements the data generate process (multivariate case with serial correlated error terms) ###
##################################################################################################################

 XY_data=function(n,Pi,beta,phi_1,mu,sigma,psi_1,garch_1,garch_2){

	uu <- mvrnorm(n=n,mu=mu,Sigma=sigma,empirical = T) # Generating multivariate normal random variables
 	u1<-uu[,1]
 	u2<-uu[,2]
 	u3<-uu[,3]
 	u4<-uu[,4]

# Generate U1 from GARCH(1,1)
  	h1<-numeric(n)
 	h2<-numeric(n)
 	h3<-numeric(n)
 	h4<-numeric(n)
 	for (k in 2:n) {

	    h1[k]<-garch_1*h1[k-1] + garch_2*(u1[k-1])^2
	    h2[k]<-garch_1*h2[k-1] + garch_2*(u2[k-1])^2
	    h3[k]<-garch_1*h3[k-1] + garch_2*(u3[k-1])^2
	    h4[k]<-garch_1*h4[k-1] + garch_2*(u4[k-1])^2

	}
 	U1<-sqrt(h1)*u1
 	U2<-sqrt(h2)*u2
 	U3<-sqrt(h3)*u3
 	U4<-sqrt(h4)*u4
 	U5 <- uu[,5]
# End of generating U1 from GARCH(1,1)

      	e1=rep(0,n)
 	e2=rep(0,n)
 	e3=rep(0,n)
 	e4=rep(0,n)
 	e1[1]=U1[1]
 	e2[1]=U2[1]
 	e3[1]=U3[1]
 	e4[1]=U4[1]

 	for(t in 2:n){
 	      e1[t]=psi_1[1]*e1[t-1]+U1[t]
 	      e2[t]=psi_1[2]*e2[t-1]+U2[t]
 	      e3[t]=psi_1[3]*e3[t-1]+U3[t]
 	      e4[t]=psi_1[4]*e4[t-1]+U4[t]
        }

 	x1=rep(0,n)
 	x2=rep(0,n)
 	x3=rep(0,n)
 	x4=rep(0,n)
 	ut=rep(0,n)
 	x1[1]=e1[1]
 	x2[1]=e2[1]
 	x3[1]=e3[1]
 	x4[1]=e4[1]
 	ut[1]=U5[1]

 	for(t in 2:n){

 	      x1[t]=Pi[1]*x1[t-1]+e1[t]
 	      x2[t]=Pi[2]*x2[t-1]+e2[t]
 	      x3[t]=Pi[3]*x3[t-1]+e3[t]
 	      x4[t]=Pi[4]*x4[t-1]+e4[t]
 	      ut[t]=phi_1*ut[t-1]+U5[t]
         }



 	 y=rep(0,n)
 	 for(t in 2:n){
 	       y[t]=beta[1]*x1[t-1]+beta[2]*x2[t-1]+beta[3]*x3[t-1]+beta[4]*x4[t-1]+ut[t]

         }

 	 XY1=cbind(x1,x2,x3,x4,y)
                   
	return(XY1)

}






#################################################################################################################################
### The code below implements both IVX-AR and IVX-KMS and return the rejection rate when u_t is AR(5) with multivariate model ###
#################################################################################################################################

 coverage=function(n,Pi,beta,phi_1,mu,sigma,psi_1,level, rep,garch_1,garch_2){

 	library(doParallel)

 	registerDoParallel(cores=8) # Eight cores are registered and used

 	phi_1_grid=seq(phi_1-0.3,phi_1+0.3,by=0.02)

 	gridnum=31
 	species=c(1:rep)

 	fit<-foreach(j=species, .combine=rbind) %dopar% {
 	library(forecast);library(copula);source("multi_function.R")
 	print(j)

 	IVX_count1=0
 	IVX_count2=0
 	IVX_count3=0
 	IVX_count4=0
 	IVX_count5=0
 
	count1_aic=0
 	count2_aic=0
 	count3_aic=0
 	count4_aic=0
 	count5_aic=0

 	count1_bic=0
 	count2_bic=0
 	count3_bic=0
 	count4_bic=0
 	count5_bic=0

 	xy=XY_data(n,Pi,beta,phi_1,mu,sigma,psi_1,garch_1,garch_2)
 	res_xy<-summary(lm(y~x1,data=as.data.frame(xy)))$resid

 	q.aic<-length(auto.arima(res_xy,stationary=T,ic="aic",allowdrift = F, allowmean = F, max.p=20,max.q=0)$coef)
 	q.bic<-length(auto.arima(res_xy,stationary=T,ic="bic",allowdrift = F, allowmean = F, max.p=20,max.q=0)$coef)

 	IVXcoefficient=rep(0,gridnum*4)
 	dim(IVXcoefficient)=c(gridnum,4)
 	IVXpvalues=rep(0,gridnum*5)
 	dim(IVXpvalues)=c(gridnum,5)
 	RSE=rep(0,gridnum)

#########################################################################

	if (q.aic<=1) {
	   phi_1_est<-arima(res_xy,order=c(1,0,0),method="ML",include.mean=F)$coef
           phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)

           for(k in 1:gridnum){

                 x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
                 x2_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]
                 x3_tilta=xy[2:n,3]-phi_1_grid[k]*xy[1:(n-1),3]
                 x4_tilta=xy[2:n,4]-phi_1_grid[k]*xy[1:(n-1),4]
                 y_tilta=xy[2:n,5]-phi_1_grid[k]*xy[1:(n-1),5]

                 xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                 ivx_results=IVX(xy_tilta,1:4,5)
                 IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                 IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                 IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                 error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                 RSE[k]=var(error)
                 
            }

            RSE_star=RSE[1]
            s_star=1
            for(s in 2:gridnum){
                 if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                  }

            }
                
            if(IVXpvalues[s_star,1]<level){count1_aic<-1}
            if(IVXpvalues[s_star,2]<level){count2_aic<-1}
            if(IVXpvalues[s_star,3]<level){count3_aic<-1}
            if(IVXpvalues[s_star,4]<level){count4_aic<-1}
            if(IVXpvalues[s_star,5]<level){count5_aic<-1}

            ivx_results_kms=IVX(xy,1:4,5)
            IVXpvalues_kms=c(as.numeric(ivx_results_kms[3,]),as.numeric(ivx_results_kms[5,1]))
            if(IVXpvalues_kms[1]<level){IVX_count1<-1}
            if(IVXpvalues_kms[2]<level){IVX_count2<-1}
            if(IVXpvalues_kms[3]<level){IVX_count3<-1}
            if(IVXpvalues_kms[4]<level){IVX_count4<-1}
            if(IVXpvalues_kms[5]<level){IVX_count5<-1}

	}

#########################################################################

	if (q.bic<=1) {

		 phi_1_est<-arima(res_xy,order=c(1,0,0),method="ML",include.mean=F)$coef
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
                       x2_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]
                       x3_tilta=xy[2:n,3]-phi_1_grid[k]*xy[1:(n-1),3]
                       x4_tilta=xy[2:n,4]-phi_1_grid[k]*xy[1:(n-1),4]
                       y_tilta=xy[2:n,5]-phi_1_grid[k]*xy[1:(n-1),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_bic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_bic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_bic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_bic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_bic<-1}

	}

#########################################################################

	if (q.aic==2) {

		 phi_1_est<-arima(res_xy,order=c(2,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		 phi_2_est<-arima(res_xy,order=c(2,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[3:n,1]-phi_1_grid[k]*xy[2:(n-1),1]-phi_2_grid[k]*xy[1:(n-2),1]
                       x2_tilta=xy[3:n,2]-phi_1_grid[k]*xy[2:(n-1),2]-phi_2_grid[k]*xy[1:(n-2),2]
                       x3_tilta=xy[3:n,3]-phi_1_grid[k]*xy[2:(n-1),3]-phi_2_grid[k]*xy[1:(n-2),3]
                       x4_tilta=xy[3:n,4]-phi_1_grid[k]*xy[2:(n-1),4]-phi_2_grid[k]*xy[1:(n-2),4]
                       y_tilta=xy[3:n,5]-phi_1_grid[k]*xy[2:(n-1),5]-phi_2_grid[k]*xy[1:(n-2),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_aic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_aic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_aic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_aic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_aic<-1}

                 ivx_results_kms=IVX(xy,1:4,5)
                 IVXpvalues_kms=c(as.numeric(ivx_results_kms[3,]),as.numeric(ivx_results_kms[5,1]))
                 if(IVXpvalues_kms[1]<level){IVX_count1<-1}
                 if(IVXpvalues_kms[2]<level){IVX_count2<-1}
                 if(IVXpvalues_kms[3]<level){IVX_count3<-1}
                 if(IVXpvalues_kms[4]<level){IVX_count4<-1}
                 if(IVXpvalues_kms[5]<level){IVX_count5<-1}

	}
 
#########################################################################

	if (q.bic==2) {

		 phi_1_est<-arima(res_xy,order=c(2,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		     phi_2_est<-arima(res_xy,order=c(2,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[3:n,1]-phi_1_grid[k]*xy[2:(n-1),1]-phi_2_grid[k]*xy[1:(n-2),1]
                       x2_tilta=xy[3:n,2]-phi_1_grid[k]*xy[2:(n-1),2]-phi_2_grid[k]*xy[1:(n-2),2]
                       x3_tilta=xy[3:n,3]-phi_1_grid[k]*xy[2:(n-1),3]-phi_2_grid[k]*xy[1:(n-2),3]
                       x4_tilta=xy[3:n,4]-phi_1_grid[k]*xy[2:(n-1),4]-phi_2_grid[k]*xy[1:(n-2),4]
                       y_tilta=xy[3:n,5]-phi_1_grid[k]*xy[2:(n-1),5]-phi_2_grid[k]*xy[1:(n-2),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_bic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_bic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_bic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_bic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_bic<-1}

	}
 
#########################################################################

	if (q.aic==3) {

		 phi_1_est<-arima(res_xy,order=c(3,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		 phi_2_est<-arima(res_xy,order=c(3,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)
		 phi_3_est<-arima(res_xy,order=c(3,0,0),method="ML",include.mean=F)$coef[3]
                 phi_3_grid=seq(phi_3_est-0.3,phi_3_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[4:n,1]-phi_1_grid[k]*xy[3:(n-1),1]-phi_2_grid[k]*xy[2:(n-2),1]-phi_3_grid[k]*xy[1:(n-3),1]
                       x2_tilta=xy[4:n,2]-phi_1_grid[k]*xy[3:(n-1),2]-phi_2_grid[k]*xy[2:(n-2),2]-phi_3_grid[k]*xy[1:(n-3),2]
                       x3_tilta=xy[4:n,3]-phi_1_grid[k]*xy[3:(n-1),3]-phi_2_grid[k]*xy[2:(n-2),3]-phi_3_grid[k]*xy[1:(n-3),3]
                       x4_tilta=xy[4:n,4]-phi_1_grid[k]*xy[3:(n-1),4]-phi_2_grid[k]*xy[2:(n-2),4]-phi_3_grid[k]*xy[1:(n-3),4]
                       y_tilta=xy[4:n,5]-phi_1_grid[k]*xy[3:(n-1),5]-phi_2_grid[k]*xy[2:(n-2),5]-phi_3_grid[k]*xy[1:(n-3),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_aic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_aic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_aic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_aic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_aic<-1}

                 ivx_results_kms=IVX(xy,1:4,5)
                 IVXpvalues_kms=c(as.numeric(ivx_results_kms[3,]),as.numeric(ivx_results_kms[5,1]))
                 if(IVXpvalues_kms[1]<level){IVX_count1<-1}
                 if(IVXpvalues_kms[2]<level){IVX_count2<-1}
                 if(IVXpvalues_kms[3]<level){IVX_count3<-1}
                 if(IVXpvalues_kms[4]<level){IVX_count4<-1}
                 if(IVXpvalues_kms[5]<level){IVX_count5<-1}

	}
 
#########################################################################

	if (q.bic==3) {

		 phi_1_est<-arima(res_xy,order=c(3,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		 phi_2_est<-arima(res_xy,order=c(3,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)
		     phi_3_est<-arima(res_xy,order=c(3,0,0),method="ML",include.mean=F)$coef[3]
                 phi_3_grid=seq(phi_3_est-0.3,phi_3_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[4:n,1]-phi_1_grid[k]*xy[3:(n-1),1]-phi_2_grid[k]*xy[2:(n-2),1]-phi_3_grid[k]*xy[1:(n-3),1]
                       x2_tilta=xy[4:n,2]-phi_1_grid[k]*xy[3:(n-1),2]-phi_2_grid[k]*xy[2:(n-2),2]-phi_3_grid[k]*xy[1:(n-3),2]
                       x3_tilta=xy[4:n,3]-phi_1_grid[k]*xy[3:(n-1),3]-phi_2_grid[k]*xy[2:(n-2),3]-phi_3_grid[k]*xy[1:(n-3),3]
                       x4_tilta=xy[4:n,4]-phi_1_grid[k]*xy[3:(n-1),4]-phi_2_grid[k]*xy[2:(n-2),4]-phi_3_grid[k]*xy[1:(n-3),4]
                       y_tilta=xy[4:n,5]-phi_1_grid[k]*xy[3:(n-1),5]-phi_2_grid[k]*xy[2:(n-2),5]-phi_3_grid[k]*xy[1:(n-3),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_bic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_bic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_bic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_bic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_bic<-1}

	}

#########################################################################

	if (q.aic==4) {

		 phi_1_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		 phi_2_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)
		 phi_3_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[3]
                 phi_3_grid=seq(phi_3_est-0.3,phi_3_est+0.3,by=0.02)
		 phi_4_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[4]
                 phi_4_grid=seq(phi_4_est-0.3,phi_4_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[5:n,1]-phi_1_grid[k]*xy[4:(n-1),1]-phi_2_grid[k]*xy[3:(n-2),1]-phi_3_grid[k]*xy[2:(n-3),1]-phi_4_grid[k]*xy[1:(n-4),1]
                       x2_tilta=xy[5:n,2]-phi_1_grid[k]*xy[4:(n-1),2]-phi_2_grid[k]*xy[3:(n-2),2]-phi_3_grid[k]*xy[2:(n-3),2]-phi_4_grid[k]*xy[1:(n-4),2]
                       x3_tilta=xy[5:n,3]-phi_1_grid[k]*xy[4:(n-1),3]-phi_2_grid[k]*xy[3:(n-2),3]-phi_3_grid[k]*xy[2:(n-3),3]-phi_4_grid[k]*xy[1:(n-4),3]
                       x4_tilta=xy[5:n,4]-phi_1_grid[k]*xy[4:(n-1),4]-phi_2_grid[k]*xy[3:(n-2),4]-phi_3_grid[k]*xy[2:(n-3),4]-phi_4_grid[k]*xy[1:(n-4),4]
                       y_tilta=xy[5:n,5]-phi_1_grid[k]*xy[4:(n-1),5]-phi_2_grid[k]*xy[3:(n-2),5]-phi_3_grid[k]*xy[2:(n-3),5]-phi_4_grid[k]*xy[1:(n-4),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_aic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_aic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_aic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_aic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_aic<-1}

                 ivx_results_kms=IVX(xy,1:4,5)
                 IVXpvalues_kms=c(as.numeric(ivx_results_kms[3,]),as.numeric(ivx_results_kms[5,1]))
                 if(IVXpvalues_kms[1]<level){IVX_count1<-1}
                 if(IVXpvalues_kms[2]<level){IVX_count2<-1}
                 if(IVXpvalues_kms[3]<level){IVX_count3<-1}
                 if(IVXpvalues_kms[4]<level){IVX_count4<-1}
                 if(IVXpvalues_kms[5]<level){IVX_count5<-1}

	}

#########################################################################
 
	if (q.bic==4) {

		 phi_1_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		 phi_2_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)
		 phi_3_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[3]
                 phi_3_grid=seq(phi_3_est-0.3,phi_3_est+0.3,by=0.02)
		 phi_4_est<-arima(res_xy,order=c(4,0,0),method="ML",include.mean=F)$coef[4]
                 phi_4_grid=seq(phi_4_est-0.3,phi_4_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[5:n,1]-phi_1_grid[k]*xy[4:(n-1),1]-phi_2_grid[k]*xy[3:(n-2),1]-phi_3_grid[k]*xy[2:(n-3),1]-phi_4_grid[k]*xy[1:(n-4),1]
                       x2_tilta=xy[5:n,2]-phi_1_grid[k]*xy[4:(n-1),2]-phi_2_grid[k]*xy[3:(n-2),2]-phi_3_grid[k]*xy[2:(n-3),2]-phi_4_grid[k]*xy[1:(n-4),2]
                       x3_tilta=xy[5:n,3]-phi_1_grid[k]*xy[4:(n-1),3]-phi_2_grid[k]*xy[3:(n-2),3]-phi_3_grid[k]*xy[2:(n-3),3]-phi_4_grid[k]*xy[1:(n-4),3]
                       x4_tilta=xy[5:n,4]-phi_1_grid[k]*xy[4:(n-1),4]-phi_2_grid[k]*xy[3:(n-2),4]-phi_3_grid[k]*xy[2:(n-3),4]-phi_4_grid[k]*xy[1:(n-4),4]
                       y_tilta=xy[5:n,5]-phi_1_grid[k]*xy[4:(n-1),5]-phi_2_grid[k]*xy[3:(n-2),5]-phi_3_grid[k]*xy[2:(n-3),5]-phi_4_grid[k]*xy[1:(n-4),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_bic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_bic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_bic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_bic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_bic<-1}

	}

#########################################################################

	if (q.aic>=5) {

		 phi_1_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		 phi_2_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)
		 phi_3_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[3]
                 phi_3_grid=seq(phi_3_est-0.3,phi_3_est+0.3,by=0.02)
		 phi_4_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[4]
                 phi_4_grid=seq(phi_4_est-0.3,phi_4_est+0.3,by=0.02)
		 phi_5_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[5]
                 phi_5_grid=seq(phi_5_est-0.3,phi_5_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
                       x2_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]
                       x3_tilta=xy[6:n,3]-phi_1_grid[k]*xy[5:(n-1),3]-phi_2_grid[k]*xy[4:(n-2),3]-phi_3_grid[k]*xy[3:(n-3),3]-phi_4_grid[k]*xy[2:(n-4),3]-phi_5_grid[k]*xy[1:(n-5),3]
                       x4_tilta=xy[6:n,4]-phi_1_grid[k]*xy[5:(n-1),4]-phi_2_grid[k]*xy[4:(n-2),4]-phi_3_grid[k]*xy[3:(n-3),4]-phi_4_grid[k]*xy[2:(n-4),4]-phi_5_grid[k]*xy[1:(n-5),4]
                       y_tilta=xy[6:n,5]-phi_1_grid[k]*xy[5:(n-1),5]-phi_2_grid[k]*xy[4:(n-2),5]-phi_3_grid[k]*xy[3:(n-3),5]-phi_4_grid[k]*xy[2:(n-4),5]-phi_5_grid[k]*xy[1:(n-5),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_aic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_aic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_aic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_aic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_aic<-1}

                 ivx_results_kms=IVX(xy,1:4,5)
                 IVXpvalues_kms=c(as.numeric(ivx_results_kms[3,]),as.numeric(ivx_results_kms[5,1]))
                 if(IVXpvalues_kms[1]<level){IVX_count1<-1}
                 if(IVXpvalues_kms[2]<level){IVX_count2<-1}
                 if(IVXpvalues_kms[3]<level){IVX_count3<-1}
                 if(IVXpvalues_kms[4]<level){IVX_count4<-1}
                 if(IVXpvalues_kms[5]<level){IVX_count5<-1}

	}

#########################################################################

	if (q.bic>=5) {

		 phi_1_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[1]
                 phi_1_grid=seq(phi_1_est-0.3,phi_1_est+0.3,by=0.02)
		 phi_2_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[2]
                 phi_2_grid=seq(phi_2_est-0.3,phi_2_est+0.3,by=0.02)
		 phi_3_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[3]
                 phi_3_grid=seq(phi_3_est-0.3,phi_3_est+0.3,by=0.02)
		 phi_4_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[4]
                 phi_4_grid=seq(phi_4_est-0.3,phi_4_est+0.3,by=0.02)
		 phi_5_est<-arima(res_xy,order=c(5,0,0),method="ML",include.mean=F)$coef[5]
                 phi_5_grid=seq(phi_5_est-0.3,phi_5_est+0.3,by=0.02)

                 for(k in 1:gridnum){

                       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
                       x2_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]
                       x3_tilta=xy[6:n,3]-phi_1_grid[k]*xy[5:(n-1),3]-phi_2_grid[k]*xy[4:(n-2),3]-phi_3_grid[k]*xy[3:(n-3),3]-phi_4_grid[k]*xy[2:(n-4),3]-phi_5_grid[k]*xy[1:(n-5),3]
                       x4_tilta=xy[6:n,4]-phi_1_grid[k]*xy[5:(n-1),4]-phi_2_grid[k]*xy[4:(n-2),4]-phi_3_grid[k]*xy[3:(n-3),4]-phi_4_grid[k]*xy[2:(n-4),4]-phi_5_grid[k]*xy[1:(n-5),4]
                       y_tilta=xy[6:n,5]-phi_1_grid[k]*xy[5:(n-1),5]-phi_2_grid[k]*xy[4:(n-2),5]-phi_3_grid[k]*xy[3:(n-3),5]-phi_4_grid[k]*xy[2:(n-4),5]-phi_5_grid[k]*xy[1:(n-5),5]

                       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
                       ivx_results=IVX(xy_tilta,1:4,5)
                       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
                       IVXpvalues[k,1:4]=as.numeric(ivx_results[3,])
                       IVXpvalues[k,5]=as.numeric(ivx_results[5,1])
                       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:4])
                       RSE[k]=var(error)
                 
                 }

                 RSE_star=RSE[1]
                 s_star=1
                 for(s in 2:gridnum){
                       if(RSE[s]<RSE_star){RSE_star=RSE[s]
                                     s_star=s
                        }

                 }
                
                 if(IVXpvalues[s_star,1]<level){count1_bic<-1}
                 if(IVXpvalues[s_star,2]<level){count2_bic<-1}
                 if(IVXpvalues[s_star,3]<level){count3_bic<-1}
                 if(IVXpvalues[s_star,4]<level){count4_bic<-1}
                 if(IVXpvalues[s_star,5]<level){count5_bic<-1}

	}

#########################################################################

	return(c(count1_aic,count2_aic,count3_aic,count4_aic,count5_aic,count1_bic,count2_bic,count3_bic,count4_bic,count5_bic,IVX_count1,IVX_count2,IVX_count3,IVX_count4,IVX_count5))
	}


 	fit
}




