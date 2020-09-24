library(rmarkdown)
############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
library(dplyr)
library(tidyr)
library(ggfortify)
library(rockchalk)
library(lattice)
library(car)
library(MASS)
library(AICcmodavg)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggridges)
library(hrbrthemes)
library(lme4)
library(GLMMadaptive)
#New dataset including stranding data from 2005-2018
lifehistoryparams<-read.csv("lifehistoryparams.csv")
summary(lifehistoryparams)
head(lifehistoryparams)
str(lifehistoryparams)
#calling which columns only to use from the dataframe
cols <- c(4:10)
lhp<-lifehistoryparams[,cols]
str(lhp)

################################################################
#make another dataset with superfamily for boxplot to check non-independence way below#
##################################################################
str(lifehistoryparams)
cols2<-c(1,4)
superfamily<-lifehistoryparams[,cols2]
str(superfamily)
superfamily$Superfamily=as.factor(superfamily$Superfamily)
str(superfamily)
#converting class from character to factor using as.factor function
lhp$habitat=as.factor(lhp$habitat)
lhp$redlistPH=as.factor(lhp$redlistPH)
str(lhp)
lhp
pairs(lhp[1:7], pch = 19,
      lower.panel=NULL)
summary(lhp)
names(lhp)
#Note: column names cannot have a space!!
lhp<-lhp %>% rename(Research.effort = values)



str(lhp)

#check for missing values
colSums(is.na(lhp))
complete.cases(lhp)
lhp2<-na.omit(lhp)
complete.cases(lhp2)
str(lhp2)
table(lhp2$habitat)
table(lhp2$redlistPH)
levels(lhp2$redlistPH)

#collapse DD and notassessed into one level of a factor using combine Levels function in rockchalk package
lhp2$redlistPH<-combineLevels(lhp2$redlistPH, levs = c("DD","notassessed"), newLabel = "notevaluated")
lhp2
head(lhp2)
str(lhp2)

#check for outliers in data using cleveland dotplots
var<-c("Research.effort","ageatsex","maxlength", "gregariousness", "strandings")
var2<-c("Research.effort","ageatsex","maxlength")
var3<-c("gregariousness", "strandings")


dotplot(as.matrix(as.matrix(lhp2[,var])),
        groups=FALSE,
        strip=strip.custom(bg='white',
                           par.strip.text=list(cex=1)),
        scales=list(x=list(relation="free", draw=TRUE),
                    y=list(relation="free",draw=FALSE)),
        col=1, cex=0.5, pch=16,
        xlab=list(label="Data range", cex=1.5),
        ylab=list(label="Data order", cex=1.5))
sum(lhp2$values==0)

#Plot covariates against each other and check collinearity, do not include response variable
names(lhp2)
Coll <- c("ageatsex", "maxlength", "gregariousness", "strandings")

#Fig 4.2
panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{usr <- par("usr"); on.exit(par(usr))
par(usr = c(0, 1, 0, 1))
r1=cor(x,y,use="pairwise.complete.obs")
r <- abs(cor(x, y,use="pairwise.complete.obs"))
txt <- format(c(r1, 0.1), digits=digits)[1]
txt <- paste(prefix, txt, sep="")
if(missing(cex.cor)) { cex <- 0.6/strwidth(txt) } else {
  cex = cex.cor}
text(0.5, 0.5, txt, cex = cex * r)}
pairs(lhp2[, Coll], lower.panel = panel.cor, cex.labels = 1.3)
####Correlation between gregariousness and strandings, maxlength and strandings


#vif function is from the car package. so install car package first
vif(glm(Research.effort~ageatsex+maxlength+gregariousness+strandings,
        family=poisson,
        data=lhp2))


#plot response variable against covariates
par(mfrow=c(2,3),mar=c(5,5,1,1))
plot(y=lhp2$Research.effort, x=lhp2$ageatsex, xlab="age at sexual", ylab= "number of publications",)
plot(Research.effort~maxlength, data=lhp2, xlab="maxmimum length", ylab="number of publications")
plot(Research.effort~gregariousness, data=lhp2, xlab="group sizes", ylab="number of publications")
plot(Research.effort~strandings, data=lhp2,  xlab="number of stranded individuals", ylab="number of publications")
boxplot(values~habitat, data=lifehistoryparams, xlab="habitat type", ylab="number of publications")

#trying to put regression line on the scatterplot
ggplot(lhp2, aes(x=strandings, y=Research.effort)) + 
  geom_point()+
  geom_smooth(method=lm, color="black")+
  labs(title="publications and age at sexual maturity",
       x="age at sexual maturity", y = "publications")+
  theme_classic()  

#separating by habitat
p <- ggplot(lhp2, aes(x=gregariousness, y=Research.effort, color=habitat, shape=habitat)) + 
  geom_point()+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(x="frequency of strandings", y = "publications")
p + theme_classic()  



glm1<-glm(Research.effort~ageatsex+maxlength+gregariousness+habitat+strandings,
          data=lhp2,
          family=poisson(link=log))
summary(glm1)
overdispersion<-glm1$deviance/glm1$df.residual
overdispersion
#alternative model to overdispersion
glm2<-glm(Research.effort~ageatsex+maxlength+gregariousness*habitat+strandings,
          data=lhp2,
          family=poisson(link=log))
summary(glm2)
overdispersion2<-glm2$deviance/glm2$df.residual
overdispersion2
#still overdispersed. interaction is not a factor

#check for influential outliers
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(cooks.distance(glm2),
     xlab="observation",
     ylab="Cook's distance",
     type="h",
     ylim=c(0,5),
     cex.lab=1.5)
abline(h=1,lty=2)
#there are no influential outliers


#try other link functions
#identity
glm3<-glm(Research.effort~ageatsex+maxlength+gregariousness+habitat+strandings,
          data=lhp2,
          family=poisson(link=identity))
summary(glm3)
overdispersion3<-glm3$deviance/glm3$df.residual
overdispersion3

#still over dispersed as overdispersion value is greater than 1.2
#square-root link
glm4<-glm(Research.effort~ageatsex+maxlength+gregariousness+habitat+strandings,
          data=lhp2,
          family=poisson(link=sqrt))
summary(glm4)
overdispersion4<-glm4$deviance/glm4$df.residual
overdispersion4

#check for non-linearity of data using pearson coefficients
#6. Non-linearity

E1 <- resid(glm1, type = "pearson")

#strandings 
xyplot(E1 ~ strandings, 
       data = lhp2,
       ylab = list("Pearson residuals", cex = 1.5),
       xlab = list("frequency of strandings", cex = 1.5),
       strip = function(bg='white', ...)
         strip.default(bg='white', ...),
       scales = list(alternating = TRUE,
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel = function(x,y){
         panel.points(x,y, col = 1, pch = 16, cex = 1.0)
         panel.loess(x,y, col = 1, lwd = 3)})

#ageatsex
xyplot(E1 ~ ageatsex, 
       data = lhp2,
       ylab = list("Pearson residuals", cex = 1.5),
       xlab = list("frequency of strandings", cex = 1.5),
       strip = function(bg='white', ...)
         strip.default(bg='white', ...),
       scales = list(alternating = TRUE,
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel = function(x,y){
         panel.points(x,y, col = 1, pch = 16, cex = 1.0)
         panel.loess(x,y, col = 1, lwd = 3)})

#maxlength 
xyplot(E1 ~ maxlength, 
       data = lhp2,
       ylab = list("Pearson residuals", cex = 1.5),
       xlab = list("frequency of strandings", cex = 1.5),
       strip = function(bg='white', ...)
         strip.default(bg='white', ...),
       scales = list(alternating = TRUE,
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel = function(x,y){
         panel.points(x,y, col = 1, pch = 16, cex = 1.0)
         panel.loess(x,y, col = 1, lwd = 3)})

#gregariousness 
xyplot(E1 ~ gregariousness, 
       data = lhp2,
       ylab = list("Pearson residuals", cex = 1.5),
       xlab = list("frequency of strandings", cex = 1.5),
       strip = function(bg='white', ...)
         strip.default(bg='white', ...),
       scales = list(alternating = TRUE,
                     x = list(relation = "free"),
                     y = list(relation = "same")),
       panel = function(x,y){
         panel.points(x,y, col = 1, pch = 16, cex = 1.0)
         panel.loess(x,y, col = 1, lwd = 3)})

#no evidence of non-linear patterns

#finally, use a negative binomial family for TRUE overdipsersion
nbglm<-glm.nb(Research.effort~ageatsex+maxlength+gregariousness+habitat+strandings,
              data=lhp2)
summary(nbglm)
odsnb <- nbglm$deviance/nbglm$df.residual
odsnb

#it is still overdispersed. We need to do model fitting

#first, plot pearson residuals against fitted values
E2 <- resid(nbglm, test='pearson')
F2 <- fitted(nbglm)

par(mfrow = c(2,3), mar = c(5,5,1,1))
plot(x = F2,
     y = E2,
     cex.lab = 1.2,
     xlab = "Fitted values",
     ylab = "Pearson residuals")
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)

#Plot residuals against all covariates in the model

plot(x = lhp2$ageatsex, 
     y = E2,
     xlab = "age at sexual maturity",
     ylab = "Pearson residuals",
     cex.lab = 1.2,
     type = "p")
abline(h = 0, lty = 2)

plot(x = lhp2$maxlength, 
     y = E2,
     xlab = "maximum length",
     ylab = "Pearson residuals",
     cex.lab = 1.2,
     type = "p")
abline(h = 0, lty = 2)

plot(x = lhp2$gregariousness, 
     y = E2,
     xlab = "group size",
     ylab = "Pearson residuals",
     cex.lab = 1.2,
     type = "p")
abline(h = 0, lty = 2)

plot(x = lhp2$strandings, 
     y = E2,
     xlab = "frequency of strandings",
     ylab = "Pearson residuals",
     cex.lab = 1.2,
     type = "p")
abline(h = 0, lty = 2)

plot(x = lhp2$habitat, 
     y = E2,
     xlab = "habitat",
     ylab = "Pearson residuals",
     cex.lab = 1.2,
     type = "p")
abline(h = 0, lty = 2)

plot(x = lhp2$redlistPH, 
     y = E2,
     xlab = "red list evaluation",
     ylab = "Pearson residuals",
     cex.lab = 1.2,
     type = "p")
abline(h = 0, lty = 2)

#pearson residuals plots show no cause for concern as all fitted values of covariates are distributed
#along the horizontal axis

#now we plot cook's distance to check for influential data points
par(mfrow=c(1,1), mar=c(5,5,2,2))
plot(cooks.distance(nbglm),
     xlab = "Observation", 
     ylab = "Cook's distance",
     type = "h", 
     ylim = c(0, 1.2),
     cex.lab =  1.5)
abline(h = 1, lty = 2)
# No evidence of influential observations.


#now create models with different combinations and sets of parameters

#1st model with no age at sexual maturity
#model name noageatsex #family=poisson
noageatsex<-glm(Research.effort~maxlength+gregariousness+habitat+strandings,
                data=lhp2, family=poisson(link=log))
summary(noageatsex)
noageatsexods<-noageatsex$deviance/noageatsex$df.residual
noageatsexods #overdispersed

#2nd model without habitat. #model name wohabitat family=poisson
wohabitat<-glm(Research.effort~ageatsex+maxlength+gregariousness+strandings,
               data=lhp2, family=poisson(link=log))
summary(wohabitat)
wohabitatods<-wohabitat$deviance/wohabitat$df.residual
wohabitatods  

#3rd model using above parameters but negative binomial family to account for overdispersion
wohabitatnb<-glm.nb(Research.effort~ageatsex+maxlength+gregariousness+strandings,
                    data=lhp2)
summary(wohabitatnb)

wohabitatnbods<-wohabitatnb$deviance/wohabitatnb$df.residual
wohabitatnbods  #still overdispersed

#4th model removing both ageatsex and habitat. still use family=poisson
glm10<-glm(Research.effort~maxlength+gregariousness+strandings,
           data=lhp2, family=poisson(link=log))
summary(glm10)  
glm10ods<-glm10$deviance/glm10$df.residual
glm10ods #overdispersed

#try above model using negative binomial family
glmnb10<-glm.nb(Research.effort~maxlength+gregariousness+strandings,
                data=lhp2)
summary(glmnb10)
glmnb10ods<-glmnb10$deviance/glmnb10$df.residual
glmnb10ods

#try poisson using fewer parameters
glm11<-glm(Research.effort~gregariousness+strandings,
           data=lhp2, family=poisson(link=log))
summary(glm11)
glm11ods<-glm11$deviance/glm11$df.residual
glm11ods
#negativebinomial of above model
glm12<-glm.nb(Research.effort~gregariousness+strandings,
              data=lhp2)
summary(glm12)
glm12ods<-glm12$deviance/glm12$df.residual
glm12ods #overdispersed

var(lhp2$Research.effort)
summary(lhp2)
#mean=8, variance=22.25


#try quasipoisson family
quasiglm<-glm(Research.effort~+ageatsex+maxlength+gregariousness+strandings,
              data=lhp2, family=quasipoisson)
summary(quasiglm)
quasiods<-quasiglm$deviance/quasiglm$df.residual
quasiods


################
################
################
#model selection using the negative binomial model with all parameters
#a negative binomial accounts for overdispersed data
##################
stepAIC(wohabitatnb)
step(wohabitatnb)
#best fit model contains 2 paramaters, max length and strandings.

model<-glm.nb(Research.effort~maxlength+strandings,
              data=lhp2,
              link="log")
summary(model)

#single term deletions using the drop1 function
drop1(model, test = "Chi")
drop1(wohabitatnb, test = "Chi")

op <- par(mfrow = c(2, 2)) 
plot(model) 
par(op)
plot(wohabitatnb)

AIC(wohabitatnb, model, glmnb10,glm11, glm12)


#create all models based on stepAIC of wohabitatnb for calculating AIC weights
model1<-glm.nb(Research.effort~ageatsex+maxlength+gregariousness+strandings,
               data=lhp2)
model2<-glm.nb(Research.effort~ageatsex+maxlength+strandings,
               data=lhp2)
model3<-glm.nb(Research.effort~maxlength+strandings,
               data=lhp2)
model4<-glm.nb(Research.effort~strandings,data=lhp2)

#comparing the models and model averaging
#since the difference AIC between models is less than 10, we do model averaging


#####
#####The best fit model is listed first.
#####
models <- list(model1, model2, model3, model4)
model.names <- c('model1', 'model2','model3', 'model4')
aictab(cand.set = models, modnames = model.names)

#model averaging using MuMIn package

options(na.action = "na.fail") ##  prevent fitting models to different datasets to prevent error on removing NAs in the dataset
#above code line is required for dredge to run

fits <- dredge(model1)
options(digits = 2)
model.sel(fits)

options(na.action = "na.omit") # set back to default
#above code sets back default

nrow(fits)

summary(model.avg(fits, subset = delta <= 2))
summary(model.avg(fits))


#to check for interaction of habitat and gregariousness and effect on publication rate

int1<-glm.nb(Research.effort~gregariousness*habitat, data=lhp2)
summary(int1)
stepAIC(int1)
#no interaction whatsoever

########################################################
########################################################

#We assume now that the model is ok and we can visualise the results


#Fig for strandings nbglm

range(lhp2$strandings) #0 - 145

MyData <- expand.grid(
  ageatsex  = mean(lhp2$ageatsex),
  maxlength = mean(lhp2$maxlength),
  gregariousness   = mean(lhp2$gregariousness),
  strandings  = seq(0, 150, length = 150))

head(MyData)

X <- model.matrix(~ ageatsex + maxlength + gregariousness + strandings, data = MyData)
head(MyData)

MyData$eta  <- X %*% coef(model1)
MyData$Pred <- exp(X %*% coef(model1))

#Calculate standard errors (SE) for predicted values
MyData$SE <- sqrt(  diag(X %*% vcov(model1) %*% t(X))  )

#And using the Pred and SE values, calculate 95% confidence intervals
MyData$SeUp <- exp(MyData$eta + 1.96 * MyData$SE)
MyData$SeLo <- exp(MyData$eta - 1.96 * MyData$SE)
head(MyData)

#And plot the observed data, fitted values and 95% CIs of the mean.
p <- ggplot()
p <- p + geom_point(data = lhp2, aes(y = Research.effort, x = strandings),shape = 16,size = 1.5)
p <- p + ylab("Research effort")
p <- p + xlab(expression(paste("Frequency of strandings")))
p <- p + ylim(-1,30)
p <- p + xlim(2.6,150)
p <- p + theme(text = element_text(size=13)) 
p <- p + theme(panel.background = element_blank())
p <- p + theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))
p <- p + theme(strip.background = element_rect(fill = "white", color = "white", size = 1))
p <- p + geom_line(data = MyData,aes(x = strandings, y = Pred), colour = "black", size = 1)
p <- p + geom_ribbon(data = MyData, aes(x = strandings, ymax = SeUp, ymin = SeLo ), alpha = 0.5)
p
#GOODFIGURE!!!!!!##


##############################
#####################################################
#Figure for nbglm. for gregariousness


range(lhp2$gregariousness) #1 - 150

MyDatagreg <- expand.grid(
  ageatsex  = mean(lhp2$ageatsex),
  maxlength = mean(lhp2$maxlength),
  strandings   = mean(lhp2$strandings),
  gregariousness  = seq(1, 150, length = 150))
head(MyDatagreg)

X2 <- model.matrix(~ ageatsex + maxlength + strandings + gregariousness, data = MyDatagreg)
head(MyDatagreg)

MyDatagreg$eta  <- X2 %*% coef(model1)
MyDatagreg$Pred <- exp(X2 %*% coef(model1))

#Calculate standard errors (SE) for predicted values
MyDatagreg$SE <- sqrt(  diag(X2 %*% vcov(model1) %*% t(X2))  )

#And using the Pred and SE values, calculate 95% confidence intervals
MyDatagreg$SeUp <- exp(MyDatagreg$eta + 1.96 * MyDatagreg$SE)
MyDatagreg$SeLo <- exp(MyDatagreg$eta - 1.96 * MyDatagreg$SE)
head(MyDatagreg)

#And plot the observed data, fitted values and 95% CIs of the mean.
pgreg <- ggplot()
pgreg <- pgreg + geom_point(data =lhp2, aes(y = Research.effort, x = gregariousness),shape = 16,size = 1.5)
pgreg <- pgreg + ylab("Research effort")
pgreg <- pgreg + xlab(expression(paste("Group size")))
pgreg <- pgreg + ylim(-1,30)
pgreg <- pgreg + xlim(0,150)
pgreg <- pgreg + theme(text = element_text(size=13)) 
pgreg <- pgreg + theme(panel.background = element_blank())
pgreg <- pgreg + theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))
pgreg <- pgreg + theme(strip.background = element_rect(fill = "white", color = "white", size = 1))
pgreg <- pgreg + geom_line(data = MyDatagreg,aes(x = gregariousness, y = Pred), colour = "black", size = 1)
pgreg <- pgreg + geom_ribbon(data = MyDatagreg, aes(x = gregariousness, ymax = SeUp, ymin = SeLo ), alpha = 0.5)
pgreg

##############################
#Fig for ageatsex

range(lhp2$ageatsex) #2.5 - 15

MyDataage <- expand.grid(
  maxlength  = mean(lhp2$maxlength),
  gregariousness = mean(lhp2$gregariousness),
  strandings   = mean(lhp2$strandings),
  ageatsex  = seq(2.5, 15, length = 15))
head(MyDataage)

X3 <- model.matrix(~ maxlength + gregariousness + strandings + ageatsex, data = MyDataage)
head(MyDataage)

MyDataage$eta  <- X3 %*% coef(model1)
MyDataage$Pred <- exp(X3 %*% coef(model1))

#Calculate standard errors (SE) for predicted values
MyDataage$SE <- sqrt(  diag(X3 %*% vcov(model1) %*% t(X3))  )

#And using the Pred and SE values, calculate 95% confidence intervals
MyDataage$SeUp <- exp(MyDataage$eta + 1.96 * MyDataage$SE)
MyDataage$SeLo <- exp(MyDataage$eta - 1.96 * MyDataage$SE)
head(MyDataage)

#And plot the observed data, fitted values and 95% CIs of the mean.
page <- ggplot()
page <- page + geom_point(data = lhp2, aes(y = Research.effort, x = ageatsex),shape = 16,size = 1.5)
page <- page + ylab("Research effort")
page <- page + xlab(expression(paste("Age at sexual maturity")))
page <- page + ylim(-1,30)
page <- page + xlim(1,15)
page <- page + theme(text = element_text(size=13)) 
page <- page + theme(panel.background = element_blank())
page <- page + theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))
page <- page + theme(strip.background = element_rect(fill = "white", color = "white", size = 1))
page <- page + geom_line(data = MyDataage,aes(x = ageatsex, y = Pred), colour = "black", size = 1)
page <- page + geom_ribbon(data = MyDataage, aes(x = ageatsex, ymax = SeUp, ymin = SeLo ), alpha = 0.5)
page
p
pgreg


##############################################################################
#Fig for maximum length


range(lhp2$maxlength) #2.1 - 25

MyDatal <- expand.grid(
  ageatsex  = mean(lhp2$ageatsex),
  gregariousness = mean(lhp2$gregariousness),
  strandings   = mean(lhp2$strandings),
  maxlength  = seq(2.1, 25, length = 25))
head(MyDatal)

X4 <- model.matrix(~ ageatsex + gregariousness + strandings + maxlength, data = MyDatal)
head(MyDatal)

MyDatal$eta  <- X4 %*% coef(model1)
MyDatal$Pred <- exp(X4 %*% coef(model1))

#Calculate standard errors (SE) for predicted values
MyDatal$SE <- sqrt(  diag(X4 %*% vcov(model1) %*% t(X4))  )

#And using the Pred and SE values, calculate 95% confidence intervals
MyDatal$SeUp <- exp(MyDatal$eta + 1.96 * MyDatal$SE)
MyDatal$SeLo <- exp(MyDatal$eta - 1.96 * MyDatal$SE)
head(MyDatal)

#And plot the observed data, fitted values and 95% CIs of the mean.
pl <- ggplot()
pl <- pl + geom_point(data = lhp2, aes(y =Research.effort, x = maxlength),shape = 16,size = 1.5)
pl <- pl + ylab("Research effort")
pl <- pl + xlab(expression(paste("Maximum length")))
pl <- pl + ylim(0,30)
pl <- pl + xlim(2,25)
pl <- pl + theme(text = element_text(size=13)) 
pl <- pl + theme(panel.background = element_blank())
pl <- pl + theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))
pl <- pl + theme(strip.background = element_rect(fill = "white", color = "white", size = 1))
pl <- pl + geom_line(data = MyDatal,aes(x = maxlength, y = Pred), colour = "black", size = 1)
pl <- pl + geom_ribbon(data = MyDatal, aes(x = maxlength, ymax = SeUp, ymin = SeLo),alpha = 0.5)

pl

#combine all 4 figure into one graph using ggpubr package
theme_set(theme_pubr())

figure <- ggarrange(page, pl, pgreg, p,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
figure




##########################################################
##########################################################
#1st part of analysis done!#
########################################################
#########################################################



####################################
#make another dataset with superfamily for boxplot to check non-independence way below#
str(lifehistoryparams)
cols2<-c(1,4)
superfamily<-lifehistoryparams[,cols2]
str(superfamily)
superfamily$Superfamily=as.factor(superfamily$Superfamily)
str(superfamily)
#converting class from character to factor using as.factor function
lhp$habitat=as.factor(lhp$habitat)
lhp$redlistPH=as.factor(lhp$redlistPH)
str(lhp)
lhp
pairs(lhp[1:7], pch = 19,
      lower.panel=NULL)
summary(lhp)





#2nd part, to model if frequency of strandings is influenced by gregariounsess 
#and unknown life history traits
#calling which columns only to use from the dataframe
cols2 <- c(1, 7,9)
A2<-lifehistoryparams[,cols2]
str(A2)
A2$Superfamily=as.factor(A2$Superfamily)
str(A2)
head(A2)


##stranding frequency as the response variable

################################################################
#RED LIST ANALYSIS
##################################################################

#set up the data format. convert redlist into binary data
lifehistoryparams<-read.csv("lifehistoryparams.csv")
summary(lifehistoryparams)
head(lifehistoryparams)
str(lifehistoryparams)
#calling which columns only to use from the dataframe
cols <- c(1:10)
lhp<-lifehistoryparams[,cols]
str(lhp)
head(lhp)

#check for missing values
colSums(is.na(lhp))
complete.cases(lhp)


#collapse DD and notassessed into one level of a factor using combine Levels function in rockchalk package
lhp$redlistPH=as.factor(lhp$redlistPH)
lhp$redlistPH<-combineLevels(lhp$redlistPH, levs = c("DD","notassessed"), 
                             newLabel = "notevaluated")
lhp

lhp$redlistPH<-combineLevels(lhp$redlistPH, levs = c("CR","VU"), newLabel = "evaluated")
levels(lhp$redlistPH)
str(lhp)

redlist.glm<-glm(redlistPH~values, data=lhp, family=binomial)
summary(redlist.glm)

#scatter plot
plot(y=lhp$redlistPH, x=lhp$values, xlab="Research effort", ylab="Red List evalutation")

#Kernel density plot
ggplot(lhp, 
       aes(x = values, 
           fill = redlistPH, xlab="Number of publications")) +
  geom_density(alpha = 0.4) +
  labs(title = "Research effort by red list evaluation")

#boxplot
ggplot(lhp, 
       aes(x = redlistPH, 
           y = values)) +
  geom_boxplot() +
  labs(title = "Research effort by red list evaluation")

#violin plot with superimposed boxplot
ggplot(lhp, 
       aes(x = redlistPH, 
           y = values)) +
  geom_violin(fill = "cornflowerblue") +
  geom_boxplot(width = .2, 
               fill = "orange",
               outlier.color = "orange",
               outlier.size = 2) + 
  labs(title = "Research effort by red list evaluation")



#ridgeline plot
str(lhp)
levels(lhp$redlistPH)
table(lhp$redlistPH)
#drop Dugongidae from plot
lhpnodugong<-filter(lhp, Superfamily != "Dugongidae")

ggplot(lhpnodugong, 
       aes(x = values, 
           y = Superfamily, 
           fill = ..x..)) +
  geom_density_ridges_gradient(alpha = 0.1, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Temp. [F]", option = "B") +
  labs(title="Research effort by taxonomic group") +
  theme_ipsum() +
  theme(legend.position="none",
        panel.spacing = unit(1, "lines"),
        strip.text.x = element_text(size = 1))

head(lhp)
str(lhp)
lhp$Superfamily=as.factor(lhp$Superfamily)
levels(lhp$Superfamily)
lhp$Superfamily=as.factor(lhp$Superfamily)





###############GLMM ANALYSIS#########################


lhp$redlistPH=as.numeric(lhp$redlistPH)
str(lhp)
head(lhp$redlistPH)

#converting redlistPH into a binary factor
lhp$redlistPH[lhp$redlistPH==1]<-0
lhp$redlistPH[lhp$redlistPH==2]<-1


#try glmer using Laplace Approximation########################
lhp.glmer<-glmer(redlistPH~values+(1|Superfamily), family=binomial, data=lhp)
summary(lhp.glmer)

#trying PQL in MASS package using Penalized quasi likelihood
LHP.PQL<-glmmPQL(redlistPH~values, random=~1|Superfamily, data=lhp, family=binomial)
summary(LHP.PQL)
#because PQL (Penalized quasi-likelihood) is the least accurate technique and suffers
#even more with small datasets like this one, let's try another glmm. see below



#GLMM using adaptive gausse hermite quadrature rule using GLMMadaptive package
glmmadaptive<-mixed_model(redlistPH~values, random=~1|Superfamily, data=lhp, family=binomial())
summary(glmmadaptive)

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################
###############STATS all done!!!!!######################
