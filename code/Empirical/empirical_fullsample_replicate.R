 ##################################################################
 ## THIS CODES CONSIDER FULL-SAMPLE SPANS FROM 1975:Q1 - 2018:Q2 ##
 ################################################################## 

 # Civillian Unemployment Rate, Seasonally Adjusted (UNRATE)
 # Real Disposable Personal Income, Percent Change from Quarter One Year Ago, Seasonally Adjusted (A067RO1Q156NBEA)
 # Gross Domestic Product, Percent Change from Preceding Period, Seasonally Adjusted Annual Rate (A191RP1Q027SBEA)
 # Effective Federal Funds Rate, Percent, Not Seasonally Adjusted (FEDFUNDS)	
 # Shares of gross domestic product: Gross private domestic investment: Fixed investment: Residential, Percent, Not Seasonally Adjusted (A011RE1Q156NBEA)
 # Industrial Production Index, Index 2012=100, Seasonally Adjusted (INDPRO)
 # Consumer Price Index for All Urban Consumers: All items less shelter, Index 1982-1984=100, Seasonally Adjusted (CUSR0000SA0L2)
 # All-Transactions House Price Index for the United States, Index 1980:Q1=100, Not Seasonally Adjusted (USSTHPI)
 # 30-Year Fixed Rate Mortgage Average in the United States, Percent, Not Seasonally Adjusted (MORTGAGE30US)
 # Gross Domestic Product: Implicit Price Deflator, Index 2012=100, Seasonally Adjusted (GDPDEF)
 # Total Reserve Balances Maintained with Federal Reserve Banks, Billions of Dollars, Not Seasonally Adjusted (RESBALNS)

 library(tseries)
 library(forecast)
 library(urca)
 library(lmtest)
 library(ggplot2)
 source("FUN_univariate.R")

 ## Load HPI data
 hpi<-read.csv("USSTHPI.csv")
 t_scale<-hpi[,1]
 hpi$per<-0
 for (i in 2:nrow(hpi)) {

 	hpi$per[i]<-(hpi[i,2]-hpi[i-1,2])/hpi[i-1,2]

 }

 price<-data.frame(time=as.Date(as.character(hpi[,1]),"%Y-%m-%d"),u=hpi[,2],v=hpi[,3])

 ## Load regressors data
 une<-subset(read.csv("UNRATE.csv"),DATE %in% t_scale)
 inc<-subset(read.csv("A067RO1Q156NBEA.csv"),DATE %in% t_scale)
 gdp<-subset(read.csv("A191RP1Q027SBEA.csv"),DATE %in% t_scale)
 cpi<-subset(read.csv("CUSR0000SA0L2.csv"),DATE %in% t_scale)
 int<-subset(read.csv("FEDFUNDS.csv"),DATE %in% t_scale)
 inv<-subset(read.csv("A011RE1Q156NBEA.csv"),DATE %in% t_scale)
 ind<-subset(read.csv("INDPRO.csv"),DATE %in% t_scale)
 mog<-subset(read.csv("MORTGAGE30US.csv"),DATE %in% t_scale)
 def<-subset(read.csv("GDPDEF.csv"),DATE %in% t_scale)
 res<-subset(read.csv("RESBALNS.csv"),DATE %in% t_scale)
 res$log_res<-log(res[,2]) # Take a log transformation on RES

 ## Unit root tests on the 10 covariates
 # 1. ADF
 adf.cpi<-ur.df(cpi[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.gdp<-ur.df(gdp[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.inc<-ur.df(inc[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.ind<-ur.df(ind[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.int<-ur.df(int[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.inv<-ur.df(inv[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.une<-ur.df(une[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.mog<-ur.df(mog[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.def<-ur.df(def[,2],type="none",lags=10,selectlags="BIC")@teststat[1]
 adf.res<-ur.df(res[,3],type="none",lags=10,selectlags="BIC")@teststat[1]

 # 2. PP test
 pp.cpi<-PP.test(cpi[,2])$statistic
 pp.gdp<-PP.test(gdp[,2])$statistic
 pp.inc<-PP.test(inc[,2])$statistic
 pp.ind<-PP.test(ind[,2])$statistic
 pp.int<-PP.test(int[,2])$statistic
 pp.inv<-PP.test(inv[,2])$statistic
 pp.une<-PP.test(une[,2])$statistic
 pp.def<-PP.test(def[,2])$statistic
 pp.mog<-PP.test(mog[,2])$statistic
 pp.res<-PP.test(res[,3])$statistic

 # 3. KPSS test
 kpss.cpi<-ur.kpss(cpi[,2],type="mu",lags="short")@teststat
 kpss.gdp<-ur.kpss(gdp[,2],type="mu",lags="short")@teststat
 kpss.inc<-ur.kpss(inc[,2],type="mu",lags="short")@teststat
 kpss.ind<-ur.kpss(ind[,2],type="mu",lags="short")@teststat
 kpss.int<-ur.kpss(int[,2],type="mu",lags="short")@teststat
 kpss.inv<-ur.kpss(inv[,2],type="mu",lags="short")@teststat
 kpss.une<-ur.kpss(une[,2],type="mu",lags="short")@teststat
 kpss.def<-ur.kpss(def[,2],type="mu",lags="short")@teststat
 kpss.mog<-ur.kpss(mog[,2],type="mu",lags="short")@teststat
 kpss.res<-ur.kpss(res[,3],type="mu",lags="short")@teststat

 # 4. DF-GLS test
 dfgls.cpi<-ur.ers(cpi[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.gdp<-ur.ers(gdp[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.inc<-ur.ers(inc[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.ind<-ur.ers(ind[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.int<-ur.ers(int[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.inv<-ur.ers(inv[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.une<-ur.ers(une[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.def<-ur.ers(def[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.mog<-ur.ers(mog[,2],type="DF-GLS",model="trend",lag.max=10)@teststat
 dfgls.res<-ur.ers(res[,3],type="DF-GLS",model="trend",lag.max=10)@teststat

 # 5. Calculate the empirical Pi
 cpi_lag<-data.frame(y=cpi[2:nrow(cpi),2],x=cpi[1:c(nrow(cpi)-1),2])
 gdp_lag<-data.frame(y=gdp[2:nrow(gdp),2],x=gdp[1:c(nrow(gdp)-1),2])
 inc_lag<-data.frame(y=inc[2:nrow(inc),2],x=inc[1:c(nrow(inc)-1),2])
 ind_lag<-data.frame(y=ind[2:nrow(ind),2],x=ind[1:c(nrow(ind)-1),2])
 int_lag<-data.frame(y=int[2:nrow(int),2],x=int[1:c(nrow(int)-1),2])
 inv_lag<-data.frame(y=inv[2:nrow(inv),2],x=inv[1:c(nrow(inv)-1),2])
 une_lag<-data.frame(y=une[2:nrow(une),2],x=une[1:c(nrow(une)-1),2])
 def_lag<-data.frame(y=def[2:nrow(def),2],x=def[1:c(nrow(def)-1),2])
 mog_lag<-data.frame(y=mog[2:nrow(mog),2],x=mog[1:c(nrow(mog)-1),2])
 res_lag<-data.frame(y=res[2:nrow(res),3],x=res[1:c(nrow(res)-1),3])

 coef.cpi<-summary(lm(y~x,data=cpi_lag))$coef[2,1]
 coef.def<-summary(lm(y~x,data=def_lag))$coef[2,1]
 coef.gdp<-summary(lm(y~x,data=gdp_lag))$coef[2,1]
 coef.inc<-summary(lm(y~x,data=inc_lag))$coef[2,1]
 coef.ind<-summary(lm(y~x,data=ind_lag))$coef[2,1]
 coef.int<-summary(lm(y~x,data=int_lag))$coef[2,1]
 coef.inv<-summary(lm(y~x,data=inv_lag))$coef[2,1]
 coef.mog<-summary(lm(y~x,data=mog_lag))$coef[2,1]
 coef.res<-summary(lm(y~x,data=res_lag))$coef[2,1]
 coef.une<-summary(lm(y~x,data=une_lag))$coef[2,1]

############################
## Replication of Table 4 ##
############################

 Table4<-cbind(round(c(coef.cpi,coef.def,coef.gdp,coef.inc,coef.ind,coef.int,coef.inv,coef.mog,coef.res,coef.une),digits=3),
               round(c(adf.cpi,adf.def,adf.gdp,adf.inc,adf.ind,adf.int,adf.inv,adf.mog,adf.res,adf.une),digits=3),
               round(c(dfgls.cpi,dfgls.def,dfgls.gdp,dfgls.inc,dfgls.ind,dfgls.int,dfgls.inv,dfgls.mog,dfgls.res,dfgls.une),digits=3),
               round(c(pp.cpi,pp.def,pp.gdp,pp.inc,pp.ind,pp.int,pp.inv,pp.mog,pp.res,pp.une),digits=3),
               round(c(kpss.cpi,kpss.def,kpss.gdp,kpss.inc,kpss.ind,kpss.int,kpss.inv,kpss.mog,kpss.res,kpss.une),digits=3))
 rownames(Table4)<-c("CPI","DEF","GDP","INC","IND","INT","INV","MOG","RES","UNE")
 colnames(Table4)<-c("Coefficient","ADF","DF-GLS","PP","KPSS")
 Table4

 # Replication of Table 4 ends here!

 ## Serial correlation test on the error terms of regressing HPI on the lag of CPI, DEF, INT and RES: Proposed Wald test, Ljung-Box, Ljung-Pierce and Breusch–Godfrey test
 # 1. OLS part: Construct yx_lag with lag in x
 yx_lag<-cbind(cpi[1:(nrow(cpi)-1),2],def[1:(nrow(def)-1),2],gdp[1:(nrow(gdp)-1),2],inc[1:(nrow(inc)-1),2],ind[1:(nrow(ind)-1),2],int[1:(nrow(int)-1),2],
               inv[1:(nrow(inv)-1),2],mog[1:(nrow(mog)-1),2],res[1:(nrow(res)-1),3],une[1:(nrow(une)-1),2],hpi[2:nrow(hpi),3])
 yx_lag<-as.data.frame(yx_lag)
 colnames(yx_lag)<-c("cpi","def","gdp","inc","ind","int","inv","mog","res","une","hpi")

 # 2. Regress HPI on the lag of CPI, DEF, INT and RES
 hpi.ols<-lm(hpi~cpi+def+int+res,data=yx_lag)
 summary(hpi.ols)

 hpi.ols.resid<-hpi.ols$residual
 residual.date<-as.Date(as.character(t_scale[-1]),"%Y-%m-%d")
 residual.dat<-data.frame(time=residual.date,u=hpi.ols.resid)

 serialtest.res<-matrix(NA,ncol=4,nrow=5)

 # (a). The proposed Wald test
 test.ar1<-arima(hpi.ols.resid,order=c(1,0,0),method="ML",include.mean = F)
 test.ar2<-arima(hpi.ols.resid,order=c(2,0,0),method="ML",include.mean = F)
 test.ar3<-arima(hpi.ols.resid,order=c(3,0,0),method="ML",include.mean = F)
 test.ar4<-arima(hpi.ols.resid,order=c(4,0,0),method="ML",include.mean = F)
 test.ar5<-arima(hpi.ols.resid,order=c(5,0,0),method="ML",include.mean = F)

 res.ar1<-test.ar1$residual ; coef.ar1<-test.ar1$coef
 res.ar2<-test.ar2$residual ; coef.ar2<-test.ar2$coef
 res.ar3<-test.ar3$residual ; coef.ar3<-test.ar3$coef
 res.ar4<-test.ar4$residual ; coef.ar4<-test.ar4$coef
 res.ar5<-test.ar5$residual ; coef.ar5<-test.ar5$coef

 res.ar<-cbind(res.ar1,res.ar2,res.ar3,res.ar4,res.ar5)

 T<-173
 W1<-(T-1-1)*coef.ar1[1]*sum(hpi.ols.resid^2)/(T-1-1)*(1/(sum(hpi.ols.resid^2*res.ar1^2)/(T-1-1)))*sum(hpi.ols.resid^2)/(T-1-1)*coef.ar1[1]

 for (q in c(2:5)) {

	u1<-matrix(hpi.ols.resid[(q+2-1):(q+2-q)],nrow=q)
	sum_u<-u1%*%t(u1)
	sum_v<-u1%*%t(u1)*res.ar[q+2,q]

 	for (t in c(q+3):T) {

	    u<-matrix(hpi.ols.resid[(t-1):(t-q)],nrow=q)
	    sum_u<-sum_u + u%*%t(u)
	    sum_v<-sum_v + u%*%t(u)*res.ar[t,q]^2

	}

	if (q==2)
	   {W<-(matrix(coef.ar2,nrow=1)%*%sum_u%*%solve(sum_v)%*%sum_u%*%matrix(coef.ar2,ncol=1))[1,1]}
	if (q==3)
	   {W<-(matrix(coef.ar3,nrow=1)%*%sum_u%*%solve(sum_v)%*%sum_u%*%matrix(coef.ar3,ncol=1))[1,1]}
	if (q==4)
	   {W<-(matrix(coef.ar4,nrow=1)%*%sum_u%*%solve(sum_v)%*%sum_u%*%matrix(coef.ar4,ncol=1))[1,1]}
	if (q==5)
	   {W<-(matrix(coef.ar5,nrow=1)%*%sum_u%*%solve(sum_v)%*%sum_u%*%matrix(coef.ar5,ncol=1))[1,1]}
	#pval<-1-pchisq(W,q)
	print(W) #; print(pval)

 }
 serialtest.res[1,1]<-W1
 serialtest.res[2,1]<-48.64896
 serialtest.res[3,1]<-47.6957
 serialtest.res[4,1]<-126.9435
 serialtest.res[5,1]<-125.4168
 
 # (b). Ljung-Box
 serialtest.res[1,2]<-Box.test(hpi.ols.resid,lag=1,type="Ljung-Box")$statistic
 serialtest.res[2,2]<-Box.test(hpi.ols.resid,lag=2,type="Ljung-Box")$statistic
 serialtest.res[3,2]<-Box.test(hpi.ols.resid,lag=3,type="Ljung-Box")$statistic
 serialtest.res[4,2]<-Box.test(hpi.ols.resid,lag=4,type="Ljung-Box")$statistic
 serialtest.res[5,2]<-Box.test(hpi.ols.resid,lag=5,type="Ljung-Box")$statistic

 # (c). Ljung-Pierce
 serialtest.res[1,3]<-Box.test(hpi.ols.resid,lag=1,type="Box-Pierce")$statistic
 serialtest.res[2,3]<-Box.test(hpi.ols.resid,lag=2,type="Box-Pierce")$statistic
 serialtest.res[3,3]<-Box.test(hpi.ols.resid,lag=3,type="Box-Pierce")$statistic
 serialtest.res[4,3]<-Box.test(hpi.ols.resid,lag=4,type="Box-Pierce")$statistic
 serialtest.res[5,3]<-Box.test(hpi.ols.resid,lag=5,type="Box-Pierce")$statistic

 # (d). Breusch–Godfrey
 serialtest.res[1,4]<-bgtest(hpi.ols.resid[2:173]~hpi.ols.resid[1:172],order=1,type="Chisq")$statistic
 serialtest.res[2,4]<-bgtest(hpi.ols.resid[2:173]~hpi.ols.resid[1:172],order=2,type="Chisq")$statistic
 serialtest.res[3,4]<-bgtest(hpi.ols.resid[2:173]~hpi.ols.resid[1:172],order=3,type="Chisq")$statistic
 serialtest.res[4,4]<-bgtest(hpi.ols.resid[2:173]~hpi.ols.resid[1:172],order=4,type="Chisq")$statistic
 serialtest.res[5,4]<-bgtest(hpi.ols.resid[2:173]~hpi.ols.resid[1:172],order=5,type="Chisq")$statistic

 serialtest.res<-as.data.frame(round(serialtest.res,digits=3))
 colnames(serialtest.res)<-c("Wald","LjungBox","BoxPierce","BreuschGodfrey")

############################
## Replication of Table 5 ##
############################
 Table5<-serialtest.res
 Table5

 ## Replication of Table 5 ends here!

 ## Choice of order (q) by AIC, AICc and BIC

 hpi.ols.resid<-hpi.ols$residual
 residual.date<-as.Date(as.character(t_scale[-1]),"%Y-%m-%d")
 residual.dat<-data.frame(time=residual.date,u=hpi.ols.resid)
 
 model<-arima(residual.dat$u,order=c(1,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar1<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(2,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc 
 model_2<-attributes(model)
 ar2<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(3,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar3<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(4,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar4<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(5,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar5<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(6,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar6<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(7,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar7<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(8,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar8<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(9,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar9<-c(model$aic,model_2$aicc,model_2$bic)

 model<-arima(residual.dat$u,order=c(10,0,0),method="ML",include.mean = F)
 npar <- length(model$coef) + 1
 nstar <- length(model$residuals) - model$arma[6] - model$arma[7] * model$arma[5]
 bic <- model$aic + npar * (log(nstar) - 2)
 aicc <- model$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)
 attr(model,"bic") <- bic
 attr(model,"aicc") <- aicc
 model_2<-attributes(model)
 ar10<-c(model$aic,model_2$aicc,model_2$bic)


############################
## Replication of Table 6 ##
############################

 Table6<-round(rbind(ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9,ar10),digits=3)
 Table6
 
 ## Replication of Table 6 ends here!

#############################
## Replication of Figure 5 ##
#############################

 windows(10,5)
 par(mar=par()$mar+c(-5,0,-3,0))
 p<-ggplot(price,aes(x=time)) + 
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1975-01-01","%Y-%m-%d"), xmax = as.Date("1975-03-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) + 
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1980-01-01","%Y-%m-%d"), xmax = as.Date("1980-07-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) + 
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1981-07-01","%Y-%m-%d"), xmax = as.Date("1982-11-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1990-07-01","%Y-%m-%d"), xmax = as.Date("1991-03-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("2001-03-01","%Y-%m-%d"), xmax = as.Date("2001-11-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("2007-12-01","%Y-%m-%d"), xmax = as.Date("2009-06-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   geom_line(aes(y = u,linetype="HPI",colour="HPI"),lwd=1.1) + geom_line(aes(y = v*10000,linetype="HPI Growth Rate",colour="HPI Growth Rate"),lwd=1.1) + 
   scale_y_continuous(sec.axis = sec_axis(~./10000, name = "Growth Rate of HPI")) + 
   scale_x_date(date_labels = '%m-%d-%Y', breaks = price$time[c(1,30,60,90,120,150,174)]) + 
   scale_colour_manual(values = c("black", "blue")) +  
   scale_linetype_manual(values = c("HPI"=6,"HPI Growth Rate"=1)) + 
   labs(y="HPI (1980:Q1=100)", x="",colour="",linetype="") +
   theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme(legend.position = c(0.5, 0.95))
 p


#############################
## Replication of Figure 6 ##
#############################

 windows(10,5)
 par(mar=par()$mar+c(0,0,-3,0))
 autoplot(acf(hpi.ols.resid, plot = FALSE, lag=20),main="")


#############################
## Replication of Figure 7 ##
#############################

 ar4<-arima(residual.dat$u,order=c(4,0,0),method="ML",include.mean = F)
 ar4_vt<-ar4$residual

 # Plot 7(a)
 windows(10,5)
 par(mar=par()$mar+c(0,0,-3,0))
 autoplot(acf(ar4_vt, plot = FALSE, lag=20),main="")

 # Plot 7(b)
 windows(10,5)
 par(mar=par()$mar+c(0,0,-3,0))
 autoplot(acf((ar4_vt)^2, plot = FALSE, lag=20),main="")

 garch11<-garch(ar4_vt,order=c(1,1),trace=F)

 # Plot 7(c)
 windows(10,5)
 par(mar=par()$mar+c(0,0,-3,0))
 autoplot(acf(garch11$residual[-1], plot = FALSE, lag=20),main="")

 # Plot 7(d)
 windows(10,5)
 par(mar=par()$mar+c(0,0,-3,0))
 autoplot(acf((garch11$residual[-1])^2, plot = FALSE, lag=20),main="")


#############################
## Replication of Figure 8 ##
#############################

 investment<-data.frame(time=as.Date(as.character(hpi[,1]),"%Y-%m-%d"),v=inv[,2])

 windows(10,5)
 par(mar=par()$mar+c(-3,0,-3,0))
 p<-ggplot(investment,aes(x=time,y=v)) + 
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1975-01-01","%Y-%m-%d"), xmax = as.Date("1975-03-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) + 
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1980-01-01","%Y-%m-%d"), xmax = as.Date("1980-07-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) + 
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1981-07-01","%Y-%m-%d"), xmax = as.Date("1982-11-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("1990-07-01","%Y-%m-%d"), xmax = as.Date("1991-03-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("2001-03-01","%Y-%m-%d"), xmax = as.Date("2001-11-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   annotate("rect", fill = "gray", alpha = 0.9, 
        xmin = as.Date("2007-12-01","%Y-%m-%d"), xmax = as.Date("2009-06-01","%Y-%m-%d"),
        ymin = -Inf, ymax = Inf) +
   ylab("Percent") + geom_line(lwd=1.1) + xlab("") + 
   scale_x_date(date_labels = '%m-%d-%Y', breaks = investment$time[c(1,30,60,90,120,150,174)]) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
 p

 ##################################################
 ## The codes below replicate Panel 1 in Table 7 ##
 ##################################################

 ## Table 7
 # 1. OLS 

 yx_lag<-cbind(cpi[1:(nrow(cpi)-1),2],def[1:(nrow(def)-1),2],gdp[1:(nrow(gdp)-1),2],inc[1:(nrow(inc)-1),2],ind[1:(nrow(ind)-1),2],int[1:(nrow(int)-1),2],
               inv[1:(nrow(inv)-1),2],mog[1:(nrow(mog)-1),2],res[1:(nrow(res)-1),3],une[1:(nrow(une)-1),2],hpi[2:nrow(hpi),3])
 yx_lag<-as.data.frame(yx_lag)

 colnames(yx_lag)[1]<-c("cpi")
 colnames(yx_lag)[2]<-c("def")
 colnames(yx_lag)[3]<-c("gdp")
 colnames(yx_lag)[4]<-c("inc")
 colnames(yx_lag)[5]<-c("ind")
 colnames(yx_lag)[6]<-c("int")
 colnames(yx_lag)[7]<-c("inv")
 colnames(yx_lag)[8]<-c("mog")
 colnames(yx_lag)[9]<-c("res")
 colnames(yx_lag)[10]<-c("une")
 colnames(yx_lag)[11]<-c("hpi")

 hpi.ols1<-lm(hpi~cpi,data=yx_lag)
 hpi.ols2<-lm(hpi~def,data=yx_lag)
 hpi.ols3<-lm(hpi~gdp,data=yx_lag)
 hpi.ols4<-lm(hpi~inc,data=yx_lag)
 hpi.ols5<-lm(hpi~ind,data=yx_lag)
 hpi.ols6<-lm(hpi~int,data=yx_lag)
 hpi.ols7<-lm(hpi~inv,data=yx_lag)
 hpi.ols8<-lm(hpi~mog,data=yx_lag)
 hpi.ols9<-lm(hpi~res,data=yx_lag)
 hpi.ols10<-lm(hpi~une,data=yx_lag)

 # 2. IVX-KMS
 yx<-cbind(cpi[,2],def[,2],gdp[,2],inc[,2],ind[,2],int[,2],inv[,2],mog[,2],res[,3],une[,2],hpi[,3])
 yx<-as.data.frame(yx)
 yx1<-yx
 b1<-IVX(data=yx1,xxrow=1,yyrow=11)
 b2<-IVX(data=yx1,xxrow=2,yyrow=11)
 b3<-IVX(data=yx1,xxrow=3,yyrow=11)
 b4<-IVX(data=yx1,xxrow=4,yyrow=11)
 b5<-IVX(data=yx1,xxrow=5,yyrow=11)
 b6<-IVX(data=yx1,xxrow=6,yyrow=11)
 b7<-IVX(data=yx1,xxrow=7,yyrow=11)
 b8<-IVX(data=yx1,xxrow=8,yyrow=11)
 b9<-IVX(data=yx1,xxrow=9,yyrow=11)
 b10<-IVX(data=yx1,xxrow=10,yyrow=11)

 # 3. IVX-AR
 # (i). CPI: AR(5)
 cpi.ar<-auto.arima(hpi.ols1$resid,d=0,max.q=0,ic="bic")
 phi_1<-cpi.ar$coef[1]
 phi_2<-cpi.ar$coef[2]
 phi_3<-cpi.ar$coef[3]
 phi_4<-cpi.ar$coef[4]
 phi_5<-cpi.ar$coef[5]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)
 phi_5_grid<-seq(phi_5-0.3,phi_5+0.3,by=0.02)

 gridnum=31
 xy<-yx[,c(1,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
 y_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]

	xy_tilta=cbind(x1_tilta,y_tilta)
	ivx_results=IVX(xy_tilta,1,2)
 	IVXcoefficient[k]=as.numeric(ivx_results[1])
 	IVXwald[k]=noquote(ivx_results[2])
 	IVXpvalues[k]=as.numeric(ivx_results[3])
 	error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
 	RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_1<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (ii). DEF: AR(5)
 def.ar<-auto.arima(hpi.ols2$resid,d=0,max.q=0,ic="bic")
 phi_1<-def.ar$coef[1]
 phi_2<-def.ar$coef[2]
 phi_3<-def.ar$coef[3]
 phi_4<-def.ar$coef[4]
 phi_5<-def.ar$coef[5]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)
 phi_5_grid<-seq(phi_5-0.3,phi_5+0.3,by=0.02)

 gridnum=31
 xy<-yx[,c(2,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
       y_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_2<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (iii). GDP: AR(1)
 gdp.ar<-auto.arima(hpi.ols3$resid,d=0,max.q=0,ic="bic")
 phi_1<-gdp.ar$coef
 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 gridnum=31
 xy<-yx[,c(3,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
       y_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_3<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (iv). INC: AR(5)
 inc.ar<-auto.arima(hpi.ols4$resid,d=0,max.q=0,ic="bic")
 phi_1<-inc.ar$coef[1]
 phi_2<-inc.ar$coef[2]
 phi_3<-inc.ar$coef[3]
 phi_4<-inc.ar$coef[4]
 phi_5<-inc.ar$coef[5]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)
 phi_5_grid<-seq(phi_5-0.3,phi_5+0.3,by=0.02)

 gridnum=31
 xy<-yx[,c(4,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
       y_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_4<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (v). IND: AR(5)
 ind.ar<-auto.arima(hpi.ols5$resid,d=0,max.q=0,ic="bic")
 phi_1<-ind.ar$coef[1]
 phi_2<-ind.ar$coef[2]
 phi_3<-ind.ar$coef[3]
 phi_4<-ind.ar$coef[4]
 phi_5<-ind.ar$coef[5]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)
 phi_5_grid<-seq(phi_5-0.3,phi_5+0.3,by=0.02)

 gridnum=31
 xy<-yx[,c(5,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
       y_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_5<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (vi). INT: AR(5)
 int.ar<-auto.arima(hpi.ols6$resid,d=0,max.q=0,ic="bic")
 phi_1<-int.ar$coef[1]
 phi_2<-int.ar$coef[2]
 phi_3<-int.ar$coef[3]
 phi_4<-int.ar$coef[4]
 phi_5<-int.ar$coef[5]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)
 phi_5_grid<-seq(phi_5-0.3,phi_5+0.3,by=0.02)

 gridnum=31
 xy<-yx[,c(6,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
       y_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_6<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (vii). INV: AR(1)
 inv.ar<-auto.arima(hpi.ols7$resid,d=0,max.q=0,ic="bic")
 phi_1<-inv.ar$coef
 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 gridnum=31
 xy<-yx[,c(7,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
       y_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_7<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (viii). MOG: AR(5)
 mog.ar<-auto.arima(hpi.ols8$resid,d=0,max.q=0,ic="bic")
 phi_1<-mog.ar$coef[1]
 phi_2<-mog.ar$coef[2]
 phi_3<-mog.ar$coef[3]
 phi_4<-mog.ar$coef[4]
 phi_5<-mog.ar$coef[5]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)
 phi_5_grid<-seq(phi_5-0.3,phi_5+0.3,by=0.02)

 gridnum=31
 xy<-yx[,c(8,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
       y_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_8<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (ix). RES: AR(5)
 res.ar<-auto.arima(hpi.ols9$resid,d=0,max.q=0,ic="bic")
 phi_1<-res.ar$coef[1]
 phi_2<-res.ar$coef[2]
 phi_3<-res.ar$coef[3]
 phi_4<-res.ar$coef[4]
 phi_5<-res.ar$coef[5]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)
 phi_5_grid<-seq(phi_5-0.3,phi_5+0.3,by=0.02)

 gridnum=31
 xy<-yx[,c(9,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[6:n,1]-phi_1_grid[k]*xy[5:(n-1),1]-phi_2_grid[k]*xy[4:(n-2),1]-phi_3_grid[k]*xy[3:(n-3),1]-phi_4_grid[k]*xy[2:(n-4),1]-phi_5_grid[k]*xy[1:(n-5),1]
       y_tilta=xy[6:n,2]-phi_1_grid[k]*xy[5:(n-1),2]-phi_2_grid[k]*xy[4:(n-2),2]-phi_3_grid[k]*xy[3:(n-3),2]-phi_4_grid[k]*xy[2:(n-4),2]-phi_5_grid[k]*xy[1:(n-5),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_9<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (x). UNE: AR(1)
 une.ar<-auto.arima(hpi.ols10$resid,d=0,max.q=0,ic="bic")
 phi_1<-une.ar$coef
 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 gridnum=31
 xy<-yx[,c(10,11)]

 IVXcoefficient=rep(0,gridnum)
 dim(IVXcoefficient)=c(gridnum)
 IVXwald=rep(0,gridnum)
 IVXpvalues=rep(0,gridnum)
 dim(IVXpvalues)=c(gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
       y_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]

       xy_tilta=cbind(x1_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1,2)
       IVXcoefficient[k]=as.numeric(ivx_results[1])
       IVXwald[k]=noquote(ivx_results[2])
       IVXpvalues[k]=as.numeric(ivx_results[3])
       error=y_tilta-sum(IVXcoefficient[k]*xy_tilta[,1])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 ivx_ar_10<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # Wald test in equation (12)
 test.cpi<-arima(hpi.ols1$resid,order=c(5,0,0),method="ML",include.mean = F)
 test.def<-arima(hpi.ols2$resid,order=c(5,0,0),method="ML",include.mean = F)
 test.gdp<-arima(hpi.ols3$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.inc<-arima(hpi.ols4$resid,order=c(5,0,0),method="ML",include.mean = F)
 test.ind<-arima(hpi.ols5$resid,order=c(5,0,0),method="ML",include.mean = F)
 test.int<-arima(hpi.ols6$resid,order=c(5,0,0),method="ML",include.mean = F)
 test.inv<-arima(hpi.ols7$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.mog<-arima(hpi.ols8$resid,order=c(5,0,0),method="ML",include.mean = F)
 test.res<-arima(hpi.ols9$resid,order=c(5,0,0),method="ML",include.mean = F)
 test.une<-arima(hpi.ols10$resid,order=c(1,0,0),method="ML",include.mean = F)

 result1<-test.cpi$residual ; coef1<-test.cpi$coef
 result2<-test.def$residual ; coef2<-test.def$coef
 result3<-test.gdp$residual ; coef3<-test.gdp$coef
 result4<-test.inc$residual ; coef4<-test.inc$coef
 result5<-test.ind$residual ; coef5<-test.ind$coef
 result6<-test.int$residual ; coef6<-test.int$coef
 result7<-test.inv$residual ; coef7<-test.inv$coef
 result8<-test.mog$residual ; coef8<-test.mog$coef
 result9<-test.res$residual ; coef9<-test.res$coef
 result10<-test.une$residual ; coef10<-test.une$coef

 T<-173
 W3<-(T-1-1)*coef3[1]*sum(hpi.ols3$resid^2)/(T-1-1)*(1/(sum(hpi.ols3$resid^2*result3^2)/(T-1-1)))*sum(hpi.ols3$resid^2)/(T-1-1)*coef3
 W7<-(T-1-1)*coef7[1]*sum(hpi.ols7$resid^2)/(T-1-1)*(1/(sum(hpi.ols7$resid^2*result7^2)/(T-1-1)))*sum(hpi.ols7$resid^2)/(T-1-1)*coef7
 W10<-(T-1-1)*coef10[1]*sum(hpi.ols10$resid^2)/(T-1-1)*(1/(sum(hpi.ols10$resid^2*result10^2)/(T-1-1)))*sum(hpi.ols10$resid^2)/(T-1-1)*coef10

 q<-5
 u1<-matrix(hpi.ols1$resid[(q+2-1):(q+2-q)],nrow=q)
 u2<-matrix(hpi.ols2$resid[(q+2-1):(q+2-q)],nrow=q)
 u4<-matrix(hpi.ols4$resid[(q+2-1):(q+2-q)],nrow=q)
 u5<-matrix(hpi.ols5$resid[(q+2-1):(q+2-q)],nrow=q)
 u6<-matrix(hpi.ols6$resid[(q+2-1):(q+2-q)],nrow=q)
 u8<-matrix(hpi.ols8$resid[(q+2-1):(q+2-q)],nrow=q)
 u9<-matrix(hpi.ols9$resid[(q+2-1):(q+2-q)],nrow=q)

 sum_u1<-u1%*%t(u1)
 sum_u2<-u2%*%t(u2)
 sum_u4<-u4%*%t(u4)
 sum_u5<-u5%*%t(u5)
 sum_u6<-u6%*%t(u6)
 sum_u8<-u8%*%t(u8)
 sum_u9<-u9%*%t(u9)

 sum_v1<-u1%*%t(u1)*result1[q+2]
 sum_v2<-u2%*%t(u2)*result2[q+2]
 sum_v4<-u4%*%t(u4)*result4[q+2]
 sum_v5<-u5%*%t(u5)*result5[q+2]
 sum_v6<-u6%*%t(u6)*result6[q+2]
 sum_v8<-u8%*%t(u8)*result8[q+2]
 sum_v9<-u9%*%t(u9)*result9[q+2]

 for (t in c(q+3):T) {

	uu1<-matrix(hpi.ols1$resid[(t-1):(t-q)],nrow=q)
	sum_u1<-sum_u1 + uu1%*%t(uu1)
	sum_v1<-sum_v1 + uu1%*%t(uu1)*result1[t]^2

	uu2<-matrix(hpi.ols2$resid[(t-1):(t-q)],nrow=q)
	sum_u2<-sum_u2 + uu2%*%t(uu2)
	sum_v2<-sum_v2 + uu2%*%t(uu2)*result2[t]^2

	uu4<-matrix(hpi.ols4$resid[(t-1):(t-q)],nrow=q)
	sum_u4<-sum_u4 + uu4%*%t(uu4)
	sum_v4<-sum_v4 + uu4%*%t(uu4)*result4[t]^2

	uu5<-matrix(hpi.ols5$resid[(t-1):(t-q)],nrow=q)
	sum_u5<-sum_u5 + uu5%*%t(uu5)
	sum_v5<-sum_v5 + uu5%*%t(uu5)*result5[t]^2

	uu6<-matrix(hpi.ols6$resid[(t-1):(t-q)],nrow=q)
	sum_u6<-sum_u6 + uu6%*%t(uu6)
	sum_v6<-sum_v6 + uu6%*%t(uu6)*result6[t]^2

	uu8<-matrix(hpi.ols8$resid[(t-1):(t-q)],nrow=q)
	sum_u8<-sum_u8 + uu8%*%t(uu8)
	sum_v8<-sum_v8 + uu8%*%t(uu8)*result8[t]^2

	uu9<-matrix(hpi.ols9$resid[(t-1):(t-q)],nrow=q)
	sum_u9<-sum_u9 + uu9%*%t(uu9)
	sum_v9<-sum_v9 + uu9%*%t(uu9)*result9[t]^2

}

 W1<-matrix(coef1,nrow=1)%*%sum_u1%*%solve(sum_v1)%*%sum_u1%*%matrix(coef1,ncol=1)
 W2<-matrix(coef2,nrow=1)%*%sum_u2%*%solve(sum_v2)%*%sum_u2%*%matrix(coef2,ncol=1)
 W4<-matrix(coef4,nrow=1)%*%sum_u4%*%solve(sum_v4)%*%sum_u4%*%matrix(coef4,ncol=1)
 W5<-matrix(coef5,nrow=1)%*%sum_u5%*%solve(sum_v5)%*%sum_u5%*%matrix(coef5,ncol=1)
 W6<-matrix(coef6,nrow=1)%*%sum_u6%*%solve(sum_v6)%*%sum_u6%*%matrix(coef6,ncol=1)
 W8<-matrix(coef8,nrow=1)%*%sum_u8%*%solve(sum_v8)%*%sum_u8%*%matrix(coef8,ncol=1)
 W9<-matrix(coef9,nrow=1)%*%sum_u9%*%solve(sum_v9)%*%sum_u9%*%matrix(coef9,ncol=1)

 # Construct the result table
 univariate.result<-matrix(NA,nrow=10,ncol=9)
 univariate.result<-as.data.frame(univariate.result)
 options(scipen=999)
 univariate.result[,1]<-round(c(hpi.ols1$coef[2],hpi.ols2$coef[2],hpi.ols3$coef[2],hpi.ols4$coef[2],hpi.ols5$coef[2],hpi.ols6$coef[2],hpi.ols7$coef[2],hpi.ols8$coef[2],hpi.ols9$coef[2],hpi.ols10$coef[2]),digits=4)
 univariate.result[,2]<-round(c(summary(hpi.ols1)$coef[2,3],summary(hpi.ols2)$coef[2,3],summary(hpi.ols3)$coef[2,3],summary(hpi.ols4)$coef[2,3],summary(hpi.ols5)$coef[2,3],summary(hpi.ols6)$coef[2,3],summary(hpi.ols7)$coef[2,3],summary(hpi.ols8)$coef[2,3],summary(hpi.ols9)$coef[2,3],summary(hpi.ols10)$coef[2,3]),digits=3)
 univariate.result[,3]<-round(c(as.numeric(b1[1]),as.numeric(b2[1]),as.numeric(b3[1]),as.numeric(b4[1]),as.numeric(b5[1]),as.numeric(b6[1]),as.numeric(b7[1]),as.numeric(b8[1]),as.numeric(b9[1]),as.numeric(b10[1])),digits=4)
 univariate.result[,4]<-c(noquote(b1[2]),noquote(b2[2]),noquote(b3[2]),noquote(b4[2]),noquote(b5[2]),noquote(b6[2]),noquote(b7[2]),noquote(b8[2]),noquote(b9[2]),noquote(b10[2]))
 univariate.result[,5]<-c(ivx_ar_1[1],ivx_ar_2[1],ivx_ar_3[1],ivx_ar_4[1],ivx_ar_5[1],ivx_ar_6[1],ivx_ar_7[1],ivx_ar_8[1],ivx_ar_9[1],ivx_ar_10[1])
 univariate.result[,6]<-c(ivx_ar_1[2],ivx_ar_2[2],ivx_ar_3[2],ivx_ar_4[2],ivx_ar_5[2],ivx_ar_6[2],ivx_ar_7[2],ivx_ar_8[2],ivx_ar_9[2],ivx_ar_10[2])
 univariate.result[,7]<-c(-0.1685,0.0872,0.2084,0.1086,0.1329,0.1795,0.3888,-0.0146,-0.2263,-0.2659)  # this is delta calculated by matlab codes
 univariate.result[,8]<-c(length(auto.arima(hpi.ols1$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols2$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols3$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols4$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols5$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols6$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols7$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols8$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols9$resid,d=0,max.q=0,ic="bic")$coef),
                          length(auto.arima(hpi.ols10$resid,d=0,max.q=0,ic="bic")$coef)) # This returns the order of AR(p)
 univariate.result[,9]<-c(W1,W2,W3,W4,W5,W6,W7,W8,W9,W10)
 colnames(univariate.result)<-c("$hat{beta}_{OLS}$","$t_{OLS}$","$hat{beta}_{IVX-KMS}$","IVX-KMS Wald","$hat{beta}_{IVX-AR}$","IVX-AR Wald","$delta$","$q$","W_phi")
 rownames(univariate.result)<-c("$CPI$","$DEF$","$GDP$","$INC$","$IND$","$INT$","$INV$","$MOG$","$RES$","$UNE$")

#######################################
## Replication of Panel 1 in Table 7 ##
#######################################
 
 Table7_1<-univariate.result
 Table7_1


 # Replication of Panel 1 in Table 7 ends here!



 ##################################################
 ## The codes below replicate Panel 1 in Table 8 ##
 ##################################################

 source("FUN_multivariate.R") # Switch the source function to multivariate model

 # 1. Combination 1: INV+UNE
 yx<-cbind(inv[,2],une[,2],hpi[,3])
 yx_com1<-as.data.frame(yx)

 colnames(yx_com1)[1]<-c("inv")
 colnames(yx_com1)[2]<-c("une")
 colnames(yx_com1)[3]<-c("hpi")

 mult.ols1<-lm(hpi~.,data=yx_com1)

 com1.ar<-auto.arima(mult.ols1$resid,d=0,max.q=0,ic="bic")
 phi_1<-com1.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

 gridnum=31
 xy<-yx_com1

 IVXcoefficient<-matrix(NA,nrow=gridnum,ncol=2)
 IVXwald<-matrix(NA,nrow=gridnum,ncol=2)
 IVXjoint<-rep(NA,gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
        print(k)
 	x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
 	x2_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]
 	y_tilta=xy[2:n,3]-phi_1_grid[k]*xy[1:(n-1),3]

 	xy_tilta=cbind(x1_tilta,x2_tilta,y_tilta)
 	ivx_results=IVX(xy_tilta,1:2,3)
 	IVXcoefficient[k,]=as.numeric(ivx_results[1,])
 	IVXwald[k,]=noquote(ivx_results[2,])
 	IVXjoint[k]=noquote(ivx_results[4,1])
 	error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:2])
 	RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 method1_com1_b<-IVXcoefficient[s_star,] ; method1_com1_joint<-noquote(IVXjoint[s_star])

 IVX_KMS_com1<-IVX(data=yx_com1,xxrow=c(1,2),yyrow=3)
 method2_com1_b<-noquote(IVX_KMS_com1[1,]) ; method2_com1_joint<-noquote(IVX_KMS_com1[4,1])

 q1<-length(com1.ar$coef)

 T<-174
 test.ar1<-arima(mult.ols1$resid,order=c(1,0,0),method="ML",include.mean = F)
 res.ar1<-test.ar1$residual ; coef.ar1<-test.ar1$coef
 W1<-(T-1-1)*coef.ar1[1]*sum(mult.ols1$resid^2)/(T-1-1)*(1/(sum(mult.ols1$resid^2*res.ar1^2)/(T-1-1)))*sum(mult.ols1$resid^2)/(T-1-1)*coef.ar1[1]


 # 2. Combination 2: GDP+INC+INDINV+UNE
 yx<-cbind(gdp[,2],inc[,2],ind[,2],inv[,2],une[,2],hpi[,3])
 yx_com2<-as.data.frame(yx)

 colnames(yx_com2)[1]<-c("gdp")
 colnames(yx_com2)[2]<-c("inc")
 colnames(yx_com2)[3]<-c("ind")
 colnames(yx_com2)[4]<-c("inv")
 colnames(yx_com2)[5]<-c("une")
 colnames(yx_com2)[6]<-c("hpi")

 mult.ols2<-lm(hpi~.,data=yx_com2)

 com2.ar<-auto.arima(mult.ols2$resid,d=0,max.q=0,ic="bic")
 phi_1<-com2.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

 gridnum=31
 xy<-yx_com2

 IVXcoefficient<-matrix(NA,nrow=gridnum,ncol=5)
 IVXwald<-matrix(NA,nrow=gridnum,ncol=5)
 IVXjoint<-rep(NA,gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
       x2_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]
       x3_tilta=xy[2:n,3]-phi_1_grid[k]*xy[1:(n-1),3]
       x4_tilta=xy[2:n,4]-phi_1_grid[k]*xy[1:(n-1),4]
       x5_tilta=xy[2:n,5]-phi_1_grid[k]*xy[1:(n-1),5]
       y_tilta=xy[2:n,6]-phi_1_grid[k]*xy[1:(n-1),6]

       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,x5_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1:5,6)
       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
       IVXwald[k,]=noquote(ivx_results[2,])
       IVXjoint[k]=noquote(ivx_results[4,1])
       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:5])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 method1_com2_b<-IVXcoefficient[s_star,] ; method1_com2_joint<-noquote(IVXjoint[s_star])

 IVX_KMS_com2<-IVX(data=yx_com2,xxrow=c(1:5),yyrow=6)
 method2_com2_b<-noquote(IVX_KMS_com2[1,]) ; method2_com2_joint<-noquote(IVX_KMS_com2[4,1])

 q2<-length(com2.ar$coef)

 T<-174
 test.ar1<-arima(mult.ols2$resid,order=c(1,0,0),method="ML",include.mean = F)
 res.ar1<-test.ar1$residual ; coef.ar1<-test.ar1$coef
 W2<-(T-1-1)*coef.ar1[1]*sum(mult.ols2$resid^2)/(T-1-1)*(1/(sum(mult.ols2$resid^2*res.ar1^2)/(T-1-1)))*sum(mult.ols2$resid^2)/(T-1-1)*coef.ar1[1]


 # 3. Combination 3: CPI+DEF+INT+RES
 yx<-cbind(cpi[,2],def[,2],int[,2],res[,3],hpi[,3])
 yx_com3<-as.data.frame(yx)

 colnames(yx_com3)[1]<-c("cpi")
 colnames(yx_com3)[2]<-c("def")
 colnames(yx_com3)[3]<-c("int")
 colnames(yx_com3)[4]<-c("res")
 colnames(yx_com3)[5]<-c("hpi")

 mult.ols3<-lm(hpi~.,data=yx_com3)

 com3.ar<-auto.arima(mult.ols3$resid,d=0,max.q=0,ic="bic")
 phi_1<-com3.ar$coef[1]
 phi_2<-com3.ar$coef[2]
 phi_3<-com3.ar$coef[3]
 phi_4<-com3.ar$coef[4]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)
 phi_2_grid<-seq(phi_2-0.3,phi_2+0.3,by=0.02)
 phi_3_grid<-seq(phi_3-0.3,phi_3+0.3,by=0.02)
 phi_4_grid<-seq(phi_4-0.3,phi_4+0.3,by=0.02)

 gridnum=31
 xy<-yx_com3

 IVXcoefficient<-matrix(NA,nrow=gridnum,ncol=4)
 IVXwald<-matrix(NA,nrow=gridnum,ncol=4)
 IVXjoint<-rep(NA,gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[5:n,1]-phi_1_grid[k]*xy[4:(n-1),1]-phi_2_grid[k]*xy[3:(n-2),1]-phi_3_grid[k]*xy[2:(n-3),1]-phi_4_grid[k]*xy[1:(n-4),1]
       x2_tilta=xy[5:n,2]-phi_1_grid[k]*xy[4:(n-1),2]-phi_2_grid[k]*xy[3:(n-2),2]-phi_3_grid[k]*xy[2:(n-3),2]-phi_4_grid[k]*xy[1:(n-4),2]
       x3_tilta=xy[5:n,3]-phi_1_grid[k]*xy[4:(n-1),3]-phi_2_grid[k]*xy[3:(n-2),3]-phi_3_grid[k]*xy[2:(n-3),3]-phi_4_grid[k]*xy[1:(n-4),3]
       x4_tilta=xy[5:n,4]-phi_1_grid[k]*xy[4:(n-1),4]-phi_2_grid[k]*xy[3:(n-2),4]-phi_3_grid[k]*xy[2:(n-3),4]-phi_4_grid[k]*xy[1:(n-4),4]
       y_tilta=xy[5:n,5]-phi_1_grid[k]*xy[4:(n-1),5]-phi_2_grid[k]*xy[3:(n-2),5]-phi_3_grid[k]*xy[2:(n-3),5]-phi_4_grid[k]*xy[1:(n-4),5]

       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1:4,5)
       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
       IVXwald[k,]=noquote(ivx_results[2,])
       IVXjoint[k]=noquote(ivx_results[4,1])
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
 method1_com3_b<-IVXcoefficient[s_star,] ; method1_com3_joint<-noquote(IVXjoint[s_star])

 IVX_KMS_com3<-IVX(data=yx_com3,xxrow=c(1:4),yyrow=5)
 method2_com3_b<-noquote(IVX_KMS_com3[1,]) ; method2_com3_joint<-noquote(IVX_KMS_com3[4,1])

 q3<-length(com3.ar$coef)

 T<-174
 test.ar3<-arima(mult.ols3$resid,order=c(4,0,0),method="ML",include.mean = F)
 res.ar3<-test.ar3$residual ; coef.ar3<-test.ar3$coef

 u1<-matrix(mult.ols3$resid[(q3+2-1):(q3+2-q3)],nrow=q3)
 sum_u<-u1%*%t(u1)

 sum_v<-u1%*%t(u1)*res.ar3[q3+2]

 for (t in c(q3+3):(T-1)) {

	u<-matrix(mult.ols3$resid[(t-1):(t-q3)],nrow=q3)
	sum_u<-sum_u + u%*%t(u)
	sum_v<-sum_v + u%*%t(u)*res.ar3[t]^2

}

 W3<-matrix(coef.ar3,nrow=1)%*%sum_u%*%solve(sum_v)%*%sum_u%*%matrix(coef.ar3,ncol=1)


 # 4. Combination 4: CPI+INT+MOG
 yx<-cbind(cpi[,2],int[,2],mog[,2],hpi[,3])
 yx_com4<-as.data.frame(yx)

 colnames(yx_com4)[1]<-c("cpi")
 colnames(yx_com4)[2]<-c("int")
 colnames(yx_com4)[3]<-c("mog")
 colnames(yx_com4)[4]<-c("hpi")

 mult.ols4<-lm(hpi~.,data=yx_com4)

 com4.ar<-auto.arima(mult.ols4$resid,d=0,max.q=0,ic="bic")
 phi_1<-com4.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

 gridnum=31
 xy<-yx_com4

 IVXcoefficient<-matrix(NA,nrow=gridnum,ncol=3)
 IVXwald<-matrix(NA,nrow=gridnum,ncol=3)
 IVXjoint<-rep(NA,gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
 print(k)
 x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
 x2_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]
 x3_tilta=xy[2:n,3]-phi_1_grid[k]*xy[1:(n-1),3]
 y_tilta=xy[2:n,4]-phi_1_grid[k]*xy[1:(n-1),4]

 xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,y_tilta)
 ivx_results=IVX(xy_tilta,1:3,4)
 IVXcoefficient[k,]=as.numeric(ivx_results[1,])
 IVXwald[k,]=noquote(ivx_results[2,])
 IVXjoint[k]=noquote(ivx_results[4,1])
 error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:3])
 RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
        if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 method1_com4_b<-IVXcoefficient[s_star,] ; method1_com4_joint<-noquote(IVXjoint[s_star])

 IVX_KMS_com4<-IVX(data=yx_com4,xxrow=c(1:3),yyrow=4)
 method2_com4_b<-noquote(IVX_KMS_com4[1,]) ; method2_com4_joint<-noquote(IVX_KMS_com4[4,1])

 q4<-length(com4.ar$coef)

 T<-174
 test.ar1<-arima(mult.ols4$resid,order=c(1,0,0),method="ML",include.mean = F)
 res.ar1<-test.ar1$residual ; coef.ar1<-test.ar1$coef
 W4<-(T-1-1)*coef.ar1[1]*sum(mult.ols4$resid^2)/(T-1-1)*(1/(sum(mult.ols4$resid^2*res.ar1^2)/(T-1-1)))*sum(mult.ols4$resid^2)/(T-1-1)*coef.ar1[1]


 # 5. Combination 5: Kitchen Sink Combination
 yx<-cbind(cpi[,2],def[,2],gdp[,2],inc[,2],ind[,2],int[,2],inv[,2],mog[,2],res[,3],une[,2],hpi[,3])
 yx_com5<-as.data.frame(yx)

 colnames(yx_com5)[1]<-c("cpi")
 colnames(yx_com5)[2]<-c("def")
 colnames(yx_com5)[3]<-c("gdp")
 colnames(yx_com5)[4]<-c("inc")
 colnames(yx_com5)[5]<-c("ind")
 colnames(yx_com5)[6]<-c("int")
 colnames(yx_com5)[7]<-c("inv")
 colnames(yx_com5)[8]<-c("mog")
 colnames(yx_com5)[9]<-c("res")
 colnames(yx_com5)[10]<-c("une")
 colnames(yx_com5)[11]<-c("hpi")

 mult.ols5<-lm(hpi~.,data=yx_com5)

 com5.ar<-auto.arima(mult.ols5$resid,d=0,max.q=0,ic="bic")
 phi_1<-com5.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

 gridnum=31
 xy<-yx_com5

 IVXcoefficient<-matrix(NA,nrow=gridnum,ncol=10)
 IVXwald<-matrix(NA,nrow=gridnum,ncol=10)
 IVXjoint<-rep(NA,gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
       x2_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]
       x3_tilta=xy[2:n,3]-phi_1_grid[k]*xy[1:(n-1),3]
       x4_tilta=xy[2:n,4]-phi_1_grid[k]*xy[1:(n-1),4]
       x5_tilta=xy[2:n,5]-phi_1_grid[k]*xy[1:(n-1),5]
       x6_tilta=xy[2:n,6]-phi_1_grid[k]*xy[1:(n-1),6]
       x7_tilta=xy[2:n,7]-phi_1_grid[k]*xy[1:(n-1),7]
       x8_tilta=xy[2:n,8]-phi_1_grid[k]*xy[1:(n-1),8]
       x9_tilta=xy[2:n,9]-phi_1_grid[k]*xy[1:(n-1),9]
       x10_tilta=xy[2:n,10]-phi_1_grid[k]*xy[1:(n-1),10]
       y_tilta=xy[2:n,11]-phi_1_grid[k]*xy[1:(n-1),11]

       xy_tilta=cbind(x1_tilta,x2_tilta,x3_tilta,x4_tilta,x5_tilta,x6_tilta,x7_tilta,x8_tilta,x9_tilta,x10_tilta,y_tilta)
       ivx_results=IVX(xy_tilta,1:10,11)
       IVXcoefficient[k,]=as.numeric(ivx_results[1,])
       IVXwald[k,]=noquote(ivx_results[2,])
       IVXjoint[k]=noquote(ivx_results[4,1])
       error=y_tilta-sum(IVXcoefficient[k,]*xy_tilta[,1:10])
       RSE[k]=var(error)
 }
 RSE_star=RSE[1]
 s_star=1
 for(s in 2:gridnum){
       if(RSE[s]<RSE_star){RSE_star=RSE[s]
		s_star=s
 	}
 }
 method1_com5_b<-IVXcoefficient[s_star,] ; method1_com5_joint<-noquote(IVXjoint[s_star])

 IVX_KMS_com5<-IVX(data=yx_com5,xxrow=c(1:10),yyrow=11)
 method2_com5_b<-noquote(IVX_KMS_com5[1,]) ; method2_com5_joint<-noquote(IVX_KMS_com5[4,1])

 q5<-length(com5.ar$coef)

 T<-174
 test.ar1<-arima(mult.ols5$resid,order=c(1,0,0),method="ML",include.mean = F)
 res.ar1<-test.ar1$residual ; coef.ar1<-test.ar1$coef
 W5<-(T-1-1)*coef.ar1[1]*sum(mult.ols5$resid^2)/(T-1-1)*(1/(sum(mult.ols5$resid^2*res.ar1^2)/(T-1-1)))*sum(mult.ols5$resid^2)/(T-1-1)*coef.ar1[1]


 ## Summarize the results
 multivariate.result<-matrix(NA,ncol=13,nrow=10)
 colnames(multivariate.result)<-c("CPI","DEF","GDP","INC","IND","INT","INV","MOG","RES","UNE","JointWald","q","W_phi")
 rownames(multivariate.result)<-rep(c("IVX-AR","IVX-KMS"),5)
 multivariate.result<-as.data.frame(multivariate.result)
 multivariate.result[1,c(7,10,11,12,13)]<-c(method1_com1_b,method1_com1_joint,q1,W1)
 multivariate.result[3,c(3:5,7,10,11,12,13)]<-c(method1_com2_b,method1_com2_joint,q2,W2)
 multivariate.result[5,c(1,2,6,9,11,12,13)]<-c(method1_com3_b,method1_com3_joint,q3,W3)
 multivariate.result[7,c(1,6,8,11,12,13)]<-c(method1_com4_b,method1_com4_joint,q4,W4)
 multivariate.result[9,c(1:13)]<-c(method1_com5_b,method1_com5_joint,q5,W5)

 multivariate.result[2,c(7,10,11)]<-c(method2_com1_b,method2_com1_joint)
 multivariate.result[4,c(3:5,7,10,11)]<-c(method2_com2_b,method2_com2_joint)
 multivariate.result[6,c(1,2,6,9,11)]<-c(method2_com3_b,method2_com3_joint)
 multivariate.result[8,c(1,6,8,11)]<-c(method2_com4_b,method2_com4_joint)
 multivariate.result[10,c(1:11)]<-c(method2_com5_b,method2_com5_joint)

#######################################
## Replication of Panel 1 in Table 8 ##
#######################################
 
 Table8_1<-multivariate.result
 Table8_1

 ## Replication of Panel 1 in Table 8 ends here!!







