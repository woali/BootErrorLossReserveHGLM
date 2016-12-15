#Example - claim reserving using HGLM
# HGLM
library(hglm)
library(tweedie)
library(statmod)
library(tweedie)
library(StatMatch) ##dummy
set.seed(456)

triangle_upper=read.csv2(file="http://web.ue.katowice.pl/woali/triangle_upper.csv")
triangle_lower=read.csv2(file="http://web.ue.katowice.pl/woali/triangle_lower.csv")
X.upper = model.matrix(~as.factor(dev),data=triangle_upper) 
Z.upper = fact2dummy(as.factor(triangle_upper$origin))
X.lower = as.matrix(cbind(rep(1,length(triangle_lower$dev)),fact2dummy(as.factor(triangle_lower$dev)))) 
Z.lower = fact2dummy(as.factor(triangle_lower$origin))
p=1

claim <-triangle_upper$claim
##GLM
beta_glm=glm(claim~as.factor(triangle_upper$dev)+as.factor(triangle_upper$origin), 
family=tweedie(var.power=1,link.power=0))

beta=beta_glm$coefficients
phi <- summary(beta_glm)$disperison

Y.upper <- exp(model.matrix(~as.factor(triangle_upper$dev)+as.factor(triangle_upper$origin))%*%beta)
Y.lower <-exp(as.matrix(cbind(rep(1,length(triangle_lower$dev)),fact2dummy(as.factor(triangle_lower$dev)),
fact2dummy(as.factor(triangle_lower$origin))))%*%beta)

Ri_glm <- tapply(Y.lower, triangle_lower$origin, sum)
R_glm=sum(Y.u.lower)

##HGLM
beta_hglm=hglm(fixed=claim~as.factor(triangle_upper$dev),random=~1|as.factor(triangle_upper$origin), 
family=tweedie(var.power=1,link.power=0), rand.family=Gamma(log))

beta=beta_hglm$fixef
u=beta_hglm$ranef
v=log(beta_hglm$ranef)
phi <- beta_hglm$varFix

Y.u.lower <- exp(X.lower%*%beta+Z.lower%*%v[2:10])
Y.u.upper <- exp(X.upper%*%beta+Z.upper%*%v)

Ri_hglm <- tapply(Y.u.lower, triangle_lower$origin, sum)
R_hglm=sum(Y.u.lower)
#################	

###############################
#Residual Bootstrap HGLM RMSEP
n.sim=1000
resid= (claim-Y.u.upper)/sqrt(Y.u.upper^p) 

RMSEP_hglm <-QAPE_hglm <-NULL
R_hglmB <- R_hglmBB <- NULL

for (i in 1:n.sim) {
residB <- sample(resid, nrow(triangle_upper), replace = TRUE)
Y.u.upperB <- abs(residB * sqrt(Y.u.upper^p) + Y.u.upper) 

beta_hglmB <- hglm(fixed=Y.u.upperB~as.factor(dev),random=~1|as.factor(origin), 
family=tweedie(var.power=p,link.power=0), rand.family=Gamma(log), data=triangle_upper)

betaB=beta_hglmB$fixef
uB=beta_hglmB$ranef
vB=log(beta_hglmB$ranef)
phiB <- beta_hglmB$varFix

Y.u.lowerB <- exp(X.lower%*%beta+Z.lower%*%vB[2:10])
Y.u.lowerBB <- rtweedie(length(Y.u.lowerB), mu = c(exp(X.lower%*%beta+Z.lower%*%vB[2:10])), 
phi = beta_hglm$varFix, power = p)

Ri_hglmB <- rbind(Ri_hglmB, tapply(Y.u.lowerB, triangle_lower$origin, sum))
Ri_hglmBB <- rbind(Ri_hglmBB, tapply(Y.u.lowerBB, triangle_lower$origin, sum))

R_hglmB=rbind(R_hglmB, sum(Y.u.lowerB))
R_hglmBB=rbind(R_hglmBB, sum(Y.u.lowerBB))
}

RMSEP_hglm <-sqrt(sum((R_hglmBB - R_hglmB)^2)/n.sim)
QAPE_hglm <-quantile(abs(R_hglmBB - R_hglmB), probs=c(0.5, 0.75,0.9, 0.95))

RMSEPi_hglm <-NULL#vector(length=9) 
QAPEi_hglm <- NULL#vector(length=9)
for (k in 1:9){
RMSEPi_hglm <-rbind(RMSEPi_hglm, sqrt(sum((Ri_hglmB[,k] - Ri_hglmBB[,k])^2)/nsim))
QAPEi_hglm <- rbind(QAPEi_hglm, quantile(abs(Ri_hglmB[,k] - Ri_hglmBB[,k]), probs=c(0.5,0.75,0.9, 0.95)))
}

# Summary
summary <-cbind(rbind(as.matrix(Ri_hglm), R_hglm), rbind(cbind(RMSEPi_hglm,QAPEi_hglm), c(RMSEP_hglm,QAPE_hglm)))
colnames(summary) <- c('R', 'resRMSEP', 'resQ0.5', 'resQ0.75', 'resQ0.9', 'resQ0.95')
summary


############## Summary
EP <-as.data.frame(cbind(RMSEPi_hglm, QAPEi_hglm))
EP

par(mfrow=c(1,2))
plot(EP[,1],xaxt="n", type="b",ylim=c(10000,350000), lty=2,pch=0, ylab="",xlab="", main="RMSEPi_HGLM and Qi_HGLM for origin")
points(EP[,2],xaxt="n",type="b",lty=2,pch=1, ylab="",xlab="")
points(EP[,3],xaxt="n",type="b", lty=2,pch=2, ylab="",xlab="")
points(EP[,4],xaxt="n",type="b", lty=2,pch=3, ylab="",xlab="")
legend("topleft",c("RMSEPi","Qi_0.5","Qi_0.75","Qi_0.9"), lty=c(2,2,2,2),pch = c(0,1,2,3), cex=0.75, inset=.05, bty = "n")
labele <-c("i=1","i=2","i=3","i=4","i=5","i=6","i=7","i=8", "i=9")
mtext(labele,at=1:9,side=1)
plot(c(RMSEP_hglm, QAPE_hglm),xaxt="n",type="b",ylab="",xlab="", main=" Total RMSEP_HGLM and Q_HGLM")
labele = c("RMSEP","Q_0.5","Q_0.75","Q_0.9")
mtext(labele,at=1:4,side=1)

