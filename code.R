rm(list=ls())
set.seed(seed=2023)
## Calling packages-----------------
library(plotly)
library(lattice)
library(MASS)
library(ggplot2)
library(latticeExtra)
library(gridExtra)
library(car)
library(tactile)
library(EnvStats)
library(dplyr)
library(reshape)
library(RVAideMemoire)
library(heplots)
library(MVN)
library(biotools)
library(forecast)
library(MVQuickGraphs)
library(qqplotr)

## calling the dataset---------------
setwd("D:/ISI/Sem-2/Swagata Nandi")
B=as.data.frame(read.csv(file="cereal.csv"))
attach(B)

## EDA of the dataset---------------------
# 1. boxplot
eda1 = list()
for(i in 1:8)
{
  eda1[[i]] = bwplot(as.matrix(B[,-c(1,2)][i]),
                     xlab=paste(colnames(B)[i+2]))
}

do.call(grid.arrange, c(eda1, ncol = 4,nrow=2, 
                        top='Boxplot of Individual Variable') )
# 2. Summary
summary(B[,-c(1,2)])

# 3. Correlation heatmap
corr = data.matrix(cor(B[,-c(1,2)]) )
mel = melt(corr)
ggplot(mel, aes(X1,X2))+geom_tile(aes(fill=value)) +
  geom_text(aes(label = round(value, 1)))+
  scale_fill_gradient2(low='red',mid = 'white' ,high='maroon')+
  labs(title="Correlation Heatmap")

## Normality Checking-------------------
# Adding some jitter to the variables as there are many repeatations
cal.jit=B[,3]+runif(43,0,7.8)
fat.jit=B[,5]+runif(43,.1,.5)
protein.jit=B[,4]+runif(43,0,.5)
fiber.jit=B[,7]+runif(43,0,2)
sugar.jit=B[,9]+runif(43,0,4)
sodium.jit=B[,6]+runif(43,0,5) #for 0 values

# original data frame modified by adding jitter
A=data.frame(Brand.of.cereal=B[,1],Manufacturer=B[,2],
             Calories=cal.jit,Protein=protein.jit,Fat=fat.jit,
             Sodium=sodium.jit,Fiber=fiber.jit,Carbohydrates=B[,8],
             Sugar=sugar.jit,Potassium=B[,10])
A
attach(A)

# Density plot of individual variable after adding noise to the variables(not considering the grouping)
plist.j = list()
for(i in 1:8)
{
  plist.j[[i]] = densityplot(~scale(A[,i+2]),data=A,
                             xlab=paste(colnames(A)[i+2]),plot.points=F,rev=F)+
    layer(panel.curve(dnorm(x),lty=2))
}
do.call(grid.arrange, c(plist.j, ncol = 4,nrow=2, 
                        top='Density Plot of Indviidual Variable after adding jitter (Standardised)'))

# QQplot of individual variables
plist_qq.j = list()
for(i in 1:8){
  plist_qq.j[[i]]=qqmath(~ scale(A[,i+2]),data = A,
                         prepanel = prepanel.qqmathline,
                         panel = function(x, ...) {
                           panel.qqmathci(x, ...)
                           panel.qqmathline(x, ...)
                           panel.qqmath(x, ...)
                         },pch=19,xlab="Theortical Quantile of N(0,1)",
                         ylab="Observed Quantiles",
                         main=paste(colnames(A)[i+2]))
}

do.call(grid.arrange, c(plist_qq.j, ncol = 4,nrow=2, 
                        top="QQPlot of Individual Variables after adding jitter (Standardised) vs N(0,1)") )

# Shapiro test----------------------
Shap1=matrix(0,nrow=8,ncol=3)
Shap1[,1]=colnames(B)[3:10]
for(i in 1:8)
{
  Shap1[i,2]=shapiro.test(A[,i+2])$p
}
Shap1[,3]=ifelse(as.numeric(Shap1[,2])<0.01,"Reject","Accept")
colnames(Shap1)=c("Variable","p value","Decision")
Shap1


## Transformation of variables (not group wise)---------------------
# Boxcox Transformation
bc <- function(data, var, optimize = FALSE){
  var <- as.character(substitute(var))
  fmla <- reformulate("1", var)
  lmvar <- lm(formula = fmla, data = as.data.frame(data), x = TRUE, y = TRUE)
  boxcox(lmvar, optimize = optimize)
}

l=array(0)
l[1]=bc(A,Calories,optimize = T)$lambda
l[2]=bc(A,Protein,optimize = T)$lambda
l[3]=bc(A,Fat,optimize = T)$lambda
l[4]=bc(A,Sodium,optimize = T)$lambda
l[5]=bc(A,Fiber,optimize = T)$lambda
l[6]=bc(A,Carbohydrates,optimize = T)$lambda
l[7]=bc(A,Sugar,optimize = T)$lambda
l[8]=bc(A,Potassium,optimize = T)$lambda
l

trans=function(X,lambda)
{
  Y=as.data.frame(matrix(0,nrow=43,ncol=10))
  Y[,1]=X[,1]
  Y[,2]=X[,2]
  for(i in 1:8)
  {
    Y[,i+2]=(X[,i+2] ^ lambda[i] - 1) / lambda[i]
  }
  colnames(Y)=colnames(X)
  return(Y)
}
D=trans(A,l)
D

## Normality Checking-----------------------------
# qqplot
plist_tq0 = list()
for(i in 1:8){
  plist_tq0[[i]]=qqmath(~ scale(D[,i+2]),data = D,
                       prepanel = prepanel.qqmathline,
                       panel = function(x, ...) {
                         panel.qqmathci(x, ...)
                         panel.qqmathline(x, ...)
                         panel.qqmath(x, ...)
                       },pch=19,xlab="Theortical Quantile of N(0,1)",
                       ylab="Observed Quantiles",
                       main=paste(colnames(D)[i+2]))
}

do.call(grid.arrange, c(plist_tq0, ncol = 4,nrow=2, 
                        top="QQPlot of Transformed Individual Variables (Standardised) vs N(0,1)") )

# shapiro test
shap.wilk=function(X)
{
  Result=matrix(0,nrow=8,ncol=3)
  Result[,1]=colnames(X)[3:10]
  for(i in 1:8)
  {
    Result[i,2]=shapiro.test(X[,i+2])$p
    Result[i,3]=ifelse(as.numeric(Result[i,2])<0.01,"Reject","Accept")
  }
  return(Result)
}
shap.wilk(D)

D.new= D[,-c(1,2)]
D.new

## PCA on the original data frame-------------------
d=B[,-c(1,2)]
p=prcomp(d,center=T,scale.=T)
p$rotation
p$x
biplot(p,scale=0,main="Biplot of Principle Component 1 and 2")

# calculate total variance explained by each principal component
var_explained=p$sdev^2/sum(p$sdev^2)
var_explained
# checking how many pc are required to explain 95%variability
prop=((cumsum((p$sdev)^2))/(sum((p$sdev)^2)))*100
prop
#create scree plot
qplot(c(1:8),var_explained)+
  geom_line()+
  xlab("Principal Component")+
  ylab("Variance Explained")+
  ggtitle("Scree Plot")+
  ylim(0,1)

## Clubbing of two groups------------------------
data1=B[,-1]
data1$Manufacturer=as.factor(data1$Manufacturer)
data1
attach(data1)

# just an exploratory approach to see which groups should be clustered---------------
# variables: Potassium, Fiber, Protein (from PCA)

# 3d plot----------------
plot_ly(x=Potassium, y=Fiber, z=Protein, type="scatter3d",mode="markers", 
        color=as.factor(Manufacturer),
        title="3D Scatterplot of Potassium,Fiber and Protein")%>%
  layout(title = 'Scatter plot of top 3 variables that explains the variation of the dataset most',
         scene = list(xaxis=list(title='Potassium'),
                      yaxis=list(title='Fiber'),
                      zaxis=list(title='Protein')),
         legend = list(title=list(text='Types of Manufacturer')))
# add group2 and group3.

# New data frame-------------
data2=data1
data2$Manufacturer=as.factor(as.numeric(replace(data2$Manufacturer,which(data2$Manufacturer==3),2)))
data2

# Everything done now onwards considering 2 groups 
# EDA of group-1-----------
#1. boxplot---------
eda.g1 = list()
for(i in 1:8)
{
  eda.g1[[i]] = bwplot(as.matrix(data2[,-1])[1:17,i],
                     xlab=paste(colnames(B)[i+2]))
}

do.call(grid.arrange, c(eda.g1, ncol = 4,nrow=2, 
                        top='Boxplot of Individual Variable of Manufacturer=1') )
#2. Summary-----------
summary(data2[,-1][1:17,])

#3. Correlation heatmap----
corr = data.matrix(cor(data2[,-1][1:17,]))
mel = melt(corr)
ggplot(mel, aes(X1,X2))+geom_tile(aes(fill=value)) +
  geom_text(aes(label = round(value, 1)))+
  scale_fill_gradient2(low='red',mid = 'white' ,high='maroon')+
  labs(title="Correlation Heatmap of Manufacturer=1")

# EDA of group-2-----------
#1. boxplot---------
eda.g2 = list()
for(i in 1:8)
{
  eda.g2[[i]] = bwplot(as.matrix(data2[,-1])[18:43,i],
                       xlab=paste(colnames(B)[i+2]))
}

do.call(grid.arrange, c(eda.g2, ncol = 4,nrow=2, 
                        top='Boxplot of Individual Variable of Manufacturer=2') )
#2. Summary-----------
summary(data2[,-1][18:43,])

#3. Correlation heatmap----
corr = data.matrix(cor(data2[,-1][18:43,]))
mel = melt(corr)
ggplot(mel, aes(X1,X2))+geom_tile(aes(fill=value)) +
  geom_text(aes(label = round(value, 1)))+
  scale_fill_gradient2(low='red',mid = 'white' ,high='maroon')+
  labs(title="Correlation Heatmap of Manufacturer=2")


## Univariate Normal---------------------
# Adding jitter-----------------------------
data3=data2
data3$Calories=data2$Calories+runif(43,0,7.8)
data3$Protein=data2$Protein+runif(43,0,.5)
data3$Fat=data2$Fat+runif(43,.1,.5)
data3$Fiber=data2$Fiber+runif(43,0,2)
data3$Sodium=data2$Sodium+runif(43,0.001,0.005) # because of 0
data3$Sugar=data2$Sugar+runif(43,0.001,0.005) #because of 0
data3

# QQplot of individual variables----------------

plist_qq1 = list()
for(i in 1:8){
  plist_qq1[[i]]=qqmath(~ scale(data3[,i+1])|Manufacturer,data = data3,
                       prepanel = prepanel.qqmathline,
                       panel = function(x, ...) {
                         panel.qqmathci(x, ...)
                         panel.qqmathline(x, ...)
                         panel.qqmath(x, ...)
                       },pch=19,xlab="Theortical Quantile of N(0,1)",ylab="Observed Quantiles",
                       main=paste(colnames(data3)[i+1]))
}

do.call(grid.arrange, c(plist_qq1, ncol = 4,nrow=2, 
                        top="QQPlot of Individual Variables (Standardised) vs N(0,1)") )


# MVN test---------------------
mvn(data3[1:17,-1],mvnTest = "royston",univariateTest="SW")
mvn(data3[18:43,-1],mvnTest = "royston",univariateTest = "SW")

# Boxcox Transformations-------------------
l=array(0)
l[1]=boxcox(lm(data3$Calories~data3$Manufacturer),optimize = T)$lambda
l[2]=boxcox(lm(data3$Protein~data3$Manufacturer),optimize = T)$lambda
l[3]=boxcox(lm(data3$Fat~data3$Manufacturer),optimize = T)$lambda
l[4]=boxcox(lm(data3$Fiber~data3$Manufacturer),optimize = T)$lambda
l[5]=boxcox(lm(data3$Potassium~data3$Manufacturer),optimize = T)$lambda
l

# Transformations----------------
data4=data3
data4[,2]=(data3[,2] ^ l[1] - 1) / l[1]
data4[,3]=(data3[,3] ^ l[2] - 1) / l[2]
data4[,4]=(data3[,4] ^ l[3] - 1) / l[3]
data4[,6]=(data3[,6] ^ l[4] - 1) / l[4]
data4[,9]=(data3[,9] ^ l[5] - 1) / l[5]
data4

# Normality Test after transformation--------------
mvn(data4[1:17,-1],mvnTest = "royston",univariateTest = "SW")
mvn(data4[18:43,-1],mvnTest = "royston",univariateTest = "SW")


# gamma plot
d_i.sq.1=mahalanobis(data4[1:17,-1],colMeans(data4[1:17,-1]),cov(data4[1:17,-1]))
di<-"chisq"
dp<-list(df=8)
gg.1<-ggplot(data=as.data.frame(d_i.sq.1),mapping=aes(sample=rchisq(17,8)))+
  stat_qq_band(distribution=di,dparams=dp)+
  stat_qq_line(distribution=di,dparams=dp)+
  stat_qq_point(distribution=di,dparams=dp)+
  labs(x="Sample Mahalanobis distance for group-1",y="Theoritical quantiles of Chi square(n=17,df=8)")
gg.1

  d_i.sq.2=mahalanobis(data4[18:43,-1],colMeans(data4[18:43,-1]),cov(data4[18:43,-1]))
gg.2<-ggplot(data=as.data.frame(d_i.sq.2),mapping=aes(sample=rchisq(26,8)))+
  stat_qq_band(distribution=di,dparams=dp)+
  stat_qq_line(distribution=di,dparams=dp)+
  stat_qq_point(distribution=di,dparams=dp)+
  labs(x="Sample Mahalanobis distance for group-2",y="Theoritical quantiles of Chi square(n=26,df=8)")
gg.2

#QQplot after transformation---------------
plist_qq2 = list()
for(i in 1:8){
  plist_qq2[[i]]=qqmath(~ scale(data4[,i+1])|Manufacturer,data = data4,
                        prepanel = prepanel.qqmathline,
                        panel = function(x, ...) {
                          panel.qqmathci(x, ...)
                          panel.qqmathline(x, ...)
                          panel.qqmath(x, ...)
                        },pch=19,xlab="Theortical Quantile of N(0,1)",ylab="Observed Quantiles",
                        main=paste(colnames(data3)[i+1]))
}

do.call(grid.arrange, c(plist_qq2, ncol = 4,nrow=2, 
                        top="QQPlot of Individual Variables (Standardised) vs N(0,1)") )


# Box M-------------------
boxM(data4[,-1],data4$Manufacturer)

# QDA--------------------
model=qda(Manufacturer~.,data=data4,CV=T)
model

(1-(length(which(data4$Manufacturer==model$class))/43))*100  # APER
data.frame(original=data4$Manufacturer,predicted=model$class) # Prediction
confusionmatrix(data4$Manufacturer,model$class)

S1=(cov(data4[,-1][1:17,]))        # Covarinace matrix of Manufacturer-1
S2=(cov(data4[,-1][18:43,]))       # Covarinace matrix of Manufacturer-2
x1.bar=apply(data4[,-1][1:17,],2,mean)    # mean vector of Manufacturer-1
x2.bar=apply(data4[,-1][18:43,],2,mean)   # mean vector of Manufacturer-2
a=(solve(S1)-solve(S2))*(-0.5)            # Quadratic coefficient
b=(t(x1.bar)%*%solve(S1))-(t(x2.bar)%*%solve(S2))  # Linear Coefficient
S1
S2
x1.bar
x2.bar
a
b

sep=function(X)        # Separation function
{
  v=array(0)   # vector of qda values
  
  for(i in 1:43)
  {
    v[i]=t(as.vector(X[i,]))%*%a%*%as.vector(X[i,])+b%*%as.vector(X[i,])
  }
  Y=matrix(0,nrow=43,ncol=2)
  Y[,1]=data4$Manufacturer
  Y[,2]=v
  colnames(Y)=c("Manufacturer","Value")
  return(Y)
}

s=sep(as.matrix(data4[,-1]))   # values
s
p=c(17/43,26/43)  # prior probabilities
k=(0.5*log(det(S1)/det(S2)))+(0.5*((t(x1.bar)%*%solve(S1)%*%x1.bar)-(t(x2.bar)%*%solve(S2)%*%x2.bar)))+log(p[2]/p[1])

k                 # separation constant

p=ggplot(NULL, aes(x=data4$Manufacturer, y=as.vector(s[,2])))+
  geom_point()+ylab("Value")+xlab("Manufacturer")

p + coord_flip()+geom_hline(yintercept=k[,1],col="red") +
  ggtitle("Separation Region")    # Separation Region

ifelse(s[,2]>k[1,],"1","2")    # Prediction on training set

#lda-----------------------------------------
p1.=list()
for(i in 1:8)
{
   p1.[[i]]=xyplot(data4$Manufacturer~data4[,i+1],xlab=paste(colnames(data4)[i+1]),
                   ylab="Manufacturer",pch=19)
}

do.call(grid.arrange, c(p1., ncol = 4,nrow=2, 
                        top='Scatterplot of each Variable(Groupwise)') )

model1=lda(Manufacturer~.,data=data4)
model1

(1-(length(which(data4$Manufacturer==predict(model1)$class))/43))*100  # APER
confusionmatrix(data4$Manufacturer,predict(model1)$class)



## Factor Analysis-------------------------
#Factor Analysis
f1=factanal(data4[,-1],factors=4)
f1$loadings^2
c=apply(f1$loadings^2,1,sum)#...communality
lamda<-f1$loadings
f1$uniquenesses #uniquenesses
psi<-diag(f1$uniqueness)
s=f1$correlation
sigma=lamda%*%t(lamda)+psi
res=s-sigma #...residual matrix

f.none=factanal(data4[,-1],factors=4,rotation="none",scores="Bartlett")
lamda.none<-f.none$loadings
psi.none<-diag(f.none$uniqueness)
s.none=f.none$correlation
sigma.none=lamda.none%*%t(lamda.none)+psi.none
res.none=s.none-sigma.none #...residual matrix
c.none=apply(f.none$loadings^2,1,sum)#communality when there is no rotation done

#..specific variances or uniqueness when there is no rotation
f.var=factanal(data4[,-1],factors=4,rotation="varimax",scores="Bartlett")
lamda.var<-f.var$loadings
psi.var<-diag(f.var$uniqueness)
s.var=f.var$correlation
sigma.var=lamda.var%*%t(lamda.var)+psi.var
res.var=s.var-sigma.var #...residual matrix
c.var=apply(f.var$loadings^2,1,sum)#communality when there is an orthogonal rotation done
#uniqueness when there is an orthogonal rotation done
f.promax=factanal(data4[,-1],factors=4,rotation="promax",scores="Bartlett")
fun<-function(x,main,col)
{
  for(i in 1:3)
  {
    for(j in (i+1):4)
    {
      plot(x[,i],x[,j],xlab=paste('Factor-',i),ylab=paste('Factor-',j),
           xlim=c(-1,1),ylim=c(-1,1),main=main,col=col,pch=19)
      abline(h=0,v=0)
      text(x[,i], x[,j], labels=colnames(B)[3:10], pos=2)
    }
  }
}
par(mfrow=c(1,1))
par(mfrow=c(2,3))
fun(f.none$loadings,main="No Rotation",col="blue")
fun(f.var$loadings,main="Varimax",col="blue")
fun(f.promax$loadings,main="Promax",col="blue")


## Confidence Interval-------------------------
#Simultaneous confidence interval
par(mfrow=c(1,1))
x=data.frame(x1=data4[,6],x2=data4[,9])
mean=c(colMeans(x))
sigma.est=cov(x)
round(solve(sigma.est),3)
c=(((43-1)*8)/(43-8))*qf(0.01,df1=8,df2=35,lower.tail=F)
t=qt((0.01/4),df=42,lower.tail=F)
confidenceEllipse(mean,eig=eigen(sigma.est),n=43,p=2,alpha=0.01,xl=c(-0.5,2.5),yl=c(2,3))
#T^2 intervals
mu1.upper=mean[1]+sqrt(c)*sqrt(sigma.est[1,1]/43)
mu1.lower=mean[1]-sqrt(c)*sqrt(sigma.est[1,1]/43)
mu2.upper=mean[2]+sqrt(c)*sqrt(sigma.est[2,2]/43)
mu2.lower=mean[2]-sqrt(c)*sqrt(sigma.est[2,2]/43)
#Bonferroni's simultaneous confidence intervals
mu1.upper.b=mean[1]+t*sqrt(sigma.est[1,1]/43)
mu1.lower.b=mean[1]-t*sqrt(sigma.est[1,1]/43)
mu2.lower.b=mean[2]-t*sqrt(sigma.est[2,2]/43)
mu2.upper.b=mean[2]+t*sqrt(sigma.est[2,2]/43)
abline(h=mu2.upper)
abline(h=mu2.lower)
abline(v=mu1.lower,lty=3)
abline(v=mu1.upper,lty=3)
abline(h=mu2.upper.b,col="red")
abline(h=mu2.lower.b,col="red")
abline(v=mu1.lower.b,lty=3,col="red")
abline(v=mu1.upper.b,lty=3,col="red")

