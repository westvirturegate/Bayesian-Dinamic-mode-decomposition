base <- cbind(5*sin(1:50),20*cos(0.2*1:50),(1:50))
base <- apply(base,2,function(x){(x-min(x))/(max(x)-min(x))})

matplot(base, type = "l")



weight2 = matrix(runif(60*3),3,60)

H <- cbind(c(rep(0,30),rep(1,30)),rep(1,60))
Beta <- matrix(c(0,1,1,0,0,0),2,3)

weight2 <- H%*%Beta


dat <- (weight2)%*%t(base)

matplot(t(dat),type="l")

Y0 = (dat[,1:(dim(dat)[2]-1)])
Y1 = (dat[,2:dim(dat)[2]])

fit = SBDMDSy((Y0),(Y1),3)

fit = XSBDMDSy((Y0),(Y1),3,H = H)

rbind(rep(1,dim(Y0)[1]),rep(0,dim(Y0)[1]))
      
      
W = fit$W
phi = fit$Phi
lambda = fit$Lambda
sample.W = fit$sample.W
sample.Beta = fit$sample.Beta
Beta = fit$Beta

plot(Re(sample.W[,1,1]),main = "Re(sample.W[,1,1])")
plot(Re(sample.W[,1,2]),main = "Re(sample.W[,1,2])")
plot(Re(sample.W[,1,3]),main = "Re(sample.W[,1,3])")
matplot(Re(sample.W[,1,]),type = "l",main = "Re(sample.W[sample,1,])")
matplot(Im(sample.W[,1,]),type = "l",main = "Im(sample.W[sample,1,])")
matplot(Re(sample.W[,2,]),type = "l",main = "Re(sample.W[sample,2,])")
matplot(Im(sample.W[,2,]),type = "l",main = "Im(sample.W[sample,2,])")
matplot(Re(sample.W[,3,]),type = "l",main = "Re(sample.W[sample,3,])")
matplot(Im(sample.W[,3,]),type = "l",main = "Im(sample.W[sample,3,])")

# plot(Re(para$sample.W[,1,1]))
# plot(Re(para$sample.W[,1,2]))
# plot(Re(para$sample.W[,1,3]))
# matplot(Re(para$sample.W[,1,]),type = "l")
# matplot(Im(para$sample.W[,1,]),type = "l")
# matplot(Re(para$sample.W[,2,]),type = "l")
# matplot(Im(para$sample.W[,2,]),type = "l")
# matplot(Re(para$sample.W[,3,]),type = "l")
# matplot(Im(para$sample.W[,3,]),type = "l")

matplot(Re(sample.Beta[,1,]),type = "l", main="Re(sample.Beta[sample,1,])")
matplot(Im(sample.Beta[,1,]),type = "l", main="Im(sample.Beta[sample,1,])")
matplot(Re(sample.Beta[,2,]),type = "l", main="Re(sample.Beta[sample,2,])")
matplot(Im(sample.Beta[,2,]),type = "l", main="Im(sample.Beta[sample,2,])")


# para <-list(W = W, lambda = lambda, sample.W = sample.W,phi = phi)

matplot(t(Y0),type = "l",main = "data")
matplot(abs(phi%*%W), type = 'l',main = "abs(phi%*%W)")

# matplot(abs(para$phi%*%para$W), type = 'l')

plot(t(Y0),abs(phi%*%W))

# plot(t(Y0),abs(para$phi%*%para$W))

plot(Re(phi[,1]),type = "l",main = "Re(phi[,1])")
plot(Re(phi[,2]),type = "l",main = "Re(phi[,2])")
plot(Re(phi[,3]),type = "l",main = "Re(phi[,3])")
plot(Im(phi[,1]),type = "l",main = "Im(phi[,1])")
plot(Im(phi[,2]),type = "l",main = "Im(phi[,2])")
plot(Im(phi[,3]),type = "l",main = "Im(phi[,3])")

matplot(Re(phi),type = "l")
matplot(Im(phi),type = "l")
matplot(abs(phi),type = "l")
image(abs(W))
image(Re(W))
image(Im(W))

image(Re(Beta),main="Re(Beta)")
image(Im(Beta),main="Im(Beta)")
image(matrix(c(0,1,1,0,0,0),2,3),main = "trueBeta")

# matplot(Re(para$phi),type = "l")
# matplot(Im(para$phi),type = "l")
# matplot(abs(para$phi),type = "l")
# image(abs(para$W))
# image(Re(para$W))
# image(Im(para$W))

plot(phi[,1],type = "l")
plot(phi[,2],type = "l")
plot(phi[,3],type = "l")


matplot(base, type = "l")








dat =abs(cbind(cos(1:50)+10,5*sin(1:50/0.5)+10))
matplot(dat,type="l")

Y0 = t(dat[1:(dim(dat)[1]-1),])
Y1 = t(dat[2:dim(dat)[1],])


fit = SBDMDSy((Y0),(Y1),3)


W = fit$W
phi = fit$Phi
lambda = fit$Lambda
sample.W = fit$sample.W
plot(Re(sample.W[,1,1]))
plot(Re(sample.W[,1,2]))
plot(Re(sample.W[,1,3]))
