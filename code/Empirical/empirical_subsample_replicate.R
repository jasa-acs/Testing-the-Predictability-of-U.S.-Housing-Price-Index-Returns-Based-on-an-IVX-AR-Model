 #################################################################
 ## THIS CODES CONSIDER SUB-SAMPLE SPANS FROM 2000:Q1 - 2018:Q2 ##
 #################################################################

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

 # Load data
 hpi<-read.csv("USSTHPI.csv")[-c(1:99),] # The first 99 observations are deleted; only observations between 2000:Q1 - 2018:Q2 are kept.
 t_scale<-hpi[,1] # t_scale saves the date within the sub-sample period
 hpi$per<-0
 for (i in 2:nrow(hpi)) {

	hpi$per[i]<-(hpi[i,2]-hpi[i-1,2])/hpi[i-1,2]

}

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

 res$log_res<-log(res[,2])



 ##################################################
 ## The codes below replicate Panel 2 in Table 7 ##
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
 # (i). CPI: AR(1)
 cpi.ar<-auto.arima(hpi.ols1$resid,d=0,max.q=0,ic="bic")
 phi_1<-cpi.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

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
 ivx_ar_1<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (ii). DEF: AR(1)
 def.ar<-auto.arima(hpi.ols2$resid,d=0,max.q=0,ic="bic")
 phi_1<-def.ar$coef[1]
 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

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

 # (iv). INC: AR(1)
 inc.ar<-auto.arima(hpi.ols4$resid,d=0,max.q=0,ic="bic")
 phi_1<-inc.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

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
 ivx_ar_4<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (v). IND: AR(1)
 ind.ar<-auto.arima(hpi.ols5$resid,d=0,max.q=0,ic="bic")
 phi_1<-ind.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

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
 ivx_ar_5<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (vi). INT: AR(1)
 int.ar<-auto.arima(hpi.ols6$resid,d=0,max.q=0,ic="bic")
 phi_1<-int.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

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

 # (viii). MOG: AR(1)
 mog.ar<-auto.arima(hpi.ols8$resid,d=0,max.q=0,ic="bic")
 phi_1<-mog.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

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
 ivx_ar_8<-noquote(c(IVXcoefficient[s_star],IVXwald[s_star]))

 # (ix). RES: AR(1)
 res.ar<-auto.arima(hpi.ols9$resid,d=0,max.q=0,ic="bic")
 phi_1<-res.ar$coef[1]

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

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
 test.cpi<-arima(hpi.ols1$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.def<-arima(hpi.ols2$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.gdp<-arima(hpi.ols3$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.inc<-arima(hpi.ols4$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.ind<-arima(hpi.ols5$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.int<-arima(hpi.ols6$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.inv<-arima(hpi.ols7$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.mog<-arima(hpi.ols8$resid,order=c(1,0,0),method="ML",include.mean = F)
 test.res<-arima(hpi.ols9$resid,order=c(1,0,0),method="ML",include.mean = F)
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
 W1<-(T-1-1)*coef1[1]*sum(hpi.ols1$resid^2)/(T-1-1)*(1/(sum(hpi.ols1$resid^2*result1^2)/(T-1-1)))*sum(hpi.ols1$resid^2)/(T-1-1)*coef1
 W2<-(T-1-1)*coef2[1]*sum(hpi.ols2$resid^2)/(T-1-1)*(1/(sum(hpi.ols2$resid^2*result2^2)/(T-1-1)))*sum(hpi.ols2$resid^2)/(T-1-1)*coef2
 W3<-(T-1-1)*coef3[1]*sum(hpi.ols3$resid^2)/(T-1-1)*(1/(sum(hpi.ols3$resid^2*result3^2)/(T-1-1)))*sum(hpi.ols3$resid^2)/(T-1-1)*coef3
 W4<-(T-1-1)*coef4[1]*sum(hpi.ols4$resid^2)/(T-1-1)*(1/(sum(hpi.ols4$resid^2*result4^2)/(T-1-1)))*sum(hpi.ols4$resid^2)/(T-1-1)*coef4
 W5<-(T-1-1)*coef5[1]*sum(hpi.ols5$resid^2)/(T-1-1)*(1/(sum(hpi.ols5$resid^2*result5^2)/(T-1-1)))*sum(hpi.ols5$resid^2)/(T-1-1)*coef5
 W6<-(T-1-1)*coef6[1]*sum(hpi.ols6$resid^2)/(T-1-1)*(1/(sum(hpi.ols6$resid^2*result6^2)/(T-1-1)))*sum(hpi.ols6$resid^2)/(T-1-1)*coef6
 W7<-(T-1-1)*coef7[1]*sum(hpi.ols7$resid^2)/(T-1-1)*(1/(sum(hpi.ols7$resid^2*result7^2)/(T-1-1)))*sum(hpi.ols7$resid^2)/(T-1-1)*coef7
 W8<-(T-1-1)*coef8[1]*sum(hpi.ols8$resid^2)/(T-1-1)*(1/(sum(hpi.ols8$resid^2*result8^2)/(T-1-1)))*sum(hpi.ols8$resid^2)/(T-1-1)*coef8
 W9<-(T-1-1)*coef9[1]*sum(hpi.ols9$resid^2)/(T-1-1)*(1/(sum(hpi.ols9$resid^2*result9^2)/(T-1-1)))*sum(hpi.ols9$resid^2)/(T-1-1)*coef9
 W10<-(T-1-1)*coef10[1]*sum(hpi.ols10$resid^2)/(T-1-1)*(1/(sum(hpi.ols10$resid^2*result10^2)/(T-1-1)))*sum(hpi.ols10$resid^2)/(T-1-1)*coef10

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
## Replication of Panel 2 in Table 7 ##
#######################################
 
 Table7_2<-univariate.result
 Table7_2


 # Replication of Panel 2 in Table 7 ends here!




 ##################################################
 ## The codes below replicate Panel 2 in Table 8 ##
 ##################################################

 source("FUN_multivariate.R")

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

 T<-75
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

 T<-75
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

 phi_1_grid<-seq(phi_1-0.3,phi_1+0.3,by=0.02)

 gridnum=31
 xy<-yx_com3

 IVXcoefficient<-matrix(NA,nrow=gridnum,ncol=4)
 IVXwald<-matrix(NA,nrow=gridnum,ncol=4)
 IVXjoint<-rep(NA,gridnum)
 RSE=rep(0,gridnum) ; n<-nrow(xy)
 for(k in 1:gridnum){
       print(k)
       x1_tilta=xy[2:n,1]-phi_1_grid[k]*xy[1:(n-1),1]
       x2_tilta=xy[2:n,2]-phi_1_grid[k]*xy[1:(n-1),2]
       x3_tilta=xy[2:n,3]-phi_1_grid[k]*xy[1:(n-1),3]
       x4_tilta=xy[2:n,4]-phi_1_grid[k]*xy[1:(n-1),4]
       y_tilta=xy[2:n,5]-phi_1_grid[k]*xy[1:(n-1),5]

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

 T<-75
 test.ar1<-arima(mult.ols3$resid,order=c(1,0,0),method="ML",include.mean = F)
 res.ar1<-test.ar1$residual ; coef.ar1<-test.ar1$coef
 W3<-(T-1-1)*coef.ar1[1]*sum(mult.ols3$resid^2)/(T-1-1)*(1/(sum(mult.ols3$resid^2*res.ar1^2)/(T-1-1)))*sum(mult.ols3$resid^2)/(T-1-1)*coef.ar1[1]


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

 T<-75
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

 T<-75
 test.ar2<-arima(mult.ols5$resid,order=c(2,0,0),method="ML",include.mean = F)
 res.ar2<-test.ar2$residual ; coef.ar2<-test.ar2$coef

 u1<-matrix(mult.ols5$resid[(q5+2-1):(q5+2-q5)],nrow=q5)
 sum_u<-u1%*%t(u1)

 sum_v<-u1%*%t(u1)*res.ar2[q5+2]

 for (t in c(q5+3):(T-1)) {

	u<-matrix(mult.ols5$resid[(t-1):(t-q5)],nrow=q5)
	sum_u<-sum_u + u%*%t(u)
	sum_v<-sum_v + u%*%t(u)*res.ar2[t]^2

}

 W5<-matrix(coef.ar2,nrow=1)%*%sum_u%*%solve(sum_v)%*%sum_u%*%matrix(coef.ar2,ncol=1)



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
## Replication of Panel 2 in Table 8 ##
#######################################
 
 Table8_2<-multivariate.result
 Table8_2

 ## Replication of Panel 2 in Table 8 ends here!!

