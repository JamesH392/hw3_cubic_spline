require(ggplot2)
require(glmnet)

# generate data or upload data
n=10;
row=100; #rows of data
x=rnorm(row,10,9) #rnorm(row,0,1); #X's
x<-x/max(x)*2*pi;
y=x^3-x^2+rnorm(row,10,9); #cubic y


cubic<-cubic_spline(x,y,2);
cubic$summary$df[2]
cubic$RSS
cubic$plott
cubic$PredictY

cubic_optimal<-optimal_knots(x,y,20);
cubic_optimal$plott
cubic_optimal$optimal_knots

YY<-cubic_optimal$RSS
XX<-cubic_optimal$Dfs
plot(XX,YY,xlab = "Dofs",ylab = "RSS");

cubic_spline<-function(x,y,K){
  #parameters
  D=3; #degree 3
  #K=2; Knots(1~20);

  #create knots (equal space or by slection. defalut: equal space)
  knots=seq(min(x), max(x), length.out = K+2)  #create knots
  knots<-knots[-1];knots<-knots[-length(knots)] #knots are

  #create matrix M for linear reg.
  M<-matrix(0, nrow = length(x), ncol = 1+D+K)
  for (i in 1:(1+D+K)) {
    if(i==1){M[,i]=rep(1,length(x)); }
    else if (i>1 && i<=D+1){ M[,i]=x^(i-1);      }
    else{  M[,i]=  x*(ifelse(x>knots[ (i-1-D) ],1,0))^D }
  }

  # Please use only lsfit or lm to calculate the splines
  X=M[,2:(1+D+K)];
  # normalize data
  X = as.matrix(t((t(X)-apply(X,2,mean))/apply(X,2,sd)))
  fit<-lm(formula=y~X)
  summary<-summary(fit)
  predictY<-predict(fit)

  #plot spline
  Data<-data.frame(X=x,Y=y,predict=predictY)
  gg<-ggplot()+xlab("std X") + ylab("Y")+ geom_point(data=Data,aes(x=X,y=Y),col="gray")+geom_vline(xintercept = knots,colour="green", linetype = "longdash")+geom_smooth(data=Data,aes(x=X,y=predictY))

  # CV RSS (normalize)
  #we use ridge=0. we use it only to get RSS/RSE of LOOCV.
  fit2<-cv.glmnet(X,y,nfolds=length(x),alpha=0,grouped=FALSE)
  which_<-fit2$lambda.min #we don't leave out any X's.
  RSE <- fit2$glmnet.fit$dev.ratio[which(fit2$glmnet.fit$lambda == which_)]
  TSS <- sum((y-mean(y))^2)
  RSS<-(1-RSE)*TSS

  return(list(RSS=RSS,RSE=RSE,Predictors=X,knots=knots,PredictY=predictY,summary=summary,plott=gg))
}




optimal_knots<-function(x,y,n){
  # compare RSS vs numof knots for optimal sol
  optimal_knots=0;
  Y<-cubic_spline(x,y,0)$RSS
  Dfs<-cubic_spline(x,y,0)$summary$df[2]
  for (i in 1:n){
    if( cubic_spline(x,y,i-1)$RSS>cubic_spline(x,y,i)$RSS){
      optimal_knots=i;

    }
    Y<-append(Y,cubic_spline(x,y,i)$RSS);
    Dfs<-append(Dfs,cubic_spline(x,y,i)$summary$df[2]);
  }
  Data2<-data.frame(Knots=0:(length(Y)-1),RSS=Y)
  gg<-ggplot()+xlab("Knots") + ylab("RSS")+geom_point(data=Data2,aes(x=Knots,y=RSS),col="black")+geom_point(data=Data2,aes(x=Data2$Knots[which.min(Data2$RSS)],y=min(RSS)),col="red")

  return(list(optimal_knots=optimal_knots,RSS=Y,plott=gg,Dfs=Dfs))

}




