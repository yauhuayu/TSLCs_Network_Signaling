# ----------------------------------------------------------------------------------------------------------------------------------------- #
# * surv 
# * merge all datasets
# ------------------------------------------------------------------------------------------------------------------------------------------
library('survival')
library('affyQCReport')
library('limma')
library('Mfuzz')
library('fdrtool')
source("sub_func.r")
source("method.r")
source("method_surv.r")
source('proc_g_plot.r')


# -------- probe-Lung-specific ------ # 																					# --- annotation --- #
mark <-  read.delim("lng_cnvt.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(mark), " rows, and ", ncol(mark), "columns\n"))

m.sym = as.character(mark[,1][drop=T])
m.n.type = as.numeric(mark[,2][drop=T])
m.n.nb = as.numeric(mark[,3][drop=T])
m.n.foc = as.numeric(mark[,4][drop=T])
m.a.pbs = as.character(mark[,9][drop=T])
m.p2.pb = as.character(mark[,6][drop=T])
cbind(m.a.pbs,m.p2.pb)[25:35,]

i.foc = grep(1,m.n.type)
i.flow = grep(2,m.n.type)
e.sym = m.sym[c(i.foc,i.flow)]


# --------- read-in surv /  exprs ------------------------------- # 			

file = c('msk1','msk2','beer4573','shedden-dfci','shedden-hlm','shedden-ummi','bild3141')
dt.files = vector('list',7)

for (n in 1:2) {
	cat("Reading network expression data ......\t",file[n],"\n\n")
	
	inf = sprintf("lng_nwEpx_%s.txt",file[n])
	data <- read.table(inf, sep="\t",header=T,strip.white=T)
	
	dt.files[[n]]$sample = as.character(data[,1][drop=T])
	dt.files[[n]]$surv.st = as.numeric(data[,2][drop=T])
	dt.files[[n]]$surv.mns = as.numeric(data[,3][drop=T])
	dt.files[[n]]$met.st = as.numeric(data[,4][drop=T])
	dt.files[[n]]$met.mns = as.numeric(data[,5][drop=T])
	dt.files[[n]]$data = data.matrix(data[,6:ncol(data)])
}
for (n in 3:6) {
	cat("Reading network expression data ......\t",file[n],"\n\n")
	
	inf = sprintf("lng_nwEpx_%s.txt",file[n])
	data <- read.table(inf, sep="\t",header=T,strip.white=T)
	
	dt.files[[n]]$sample = as.character(data[,1][drop=T])
	dt.files[[n]]$surv.st = as.numeric(data[,2][drop=T])
	dt.files[[n]]$surv.mns = as.numeric(data[,3][drop=T])
	dt.files[[n]]$data = data.matrix(data[,4:ncol(data)])
}
# --------------------------------------------------------------- # 
n=7
	cat("Reading network expression data ......\t",file[n],"\n\n")
	
	inf = sprintf("lng_nwEpx_%s.txt",file[n])
	data <- read.table(inf, sep="\t",header=T,strip.white=T)
	
	dt.files[[n]]$sample = as.character(data[,1][drop=T])
	dt.files[[n]]$surv.st = as.numeric(data[,2][drop=T])
	dt.files[[n]]$surv.mns = as.numeric(data[,3][drop=T])
	dt.files[[n]]$data = data.matrix(data[,4:ncol(data)])


# ==================================================== #
# merge
# ---------------------------------------------------- #
# *** Surv *** #
epx.surv = e.sym	
surv.st = c()
surv.mns = c()
for (e in 1:6) {
	dd = dt.files[[e]]$data
	dd.1 = dd[,match(e.sym,colnames(dd),nomatch=0)]
	dd.2 = cbind(dd.1,rep(0,nrow(dd.1)))
	colnames(dd.2) = e.sym
	epx.surv = rbind(epx.surv,dd.2)
	surv.st = append(surv.st,dt.files[[e]]$surv.st)
	surv.mns = append(surv.mns,dt.files[[e]]$surv.mns)
}
e = 7 
dd = dt.files[[e]]$data
dd.1 = dd[,match(e.sym,colnames(dd),nomatch=0)]
epx.surv = rbind(epx.surv,dd.1)
surv.st = append(surv.st,dt.files[[e]]$surv.st)
surv.mns = append(surv.mns,dt.files[[e]]$surv.mns)
epx.surv = epx.surv[2:nrow(epx.surv),]
epx.surv = t(epx.surv)
k = dim(epx.surv)
epx.surv = as.numeric(epx.surv)
dim(epx.surv) = k
rownames(epx.surv) = e.sym		
	
	
# *********************************************************************************************************** #
epx.surv = epx.surv[1:68,]										# remove missing 'DDX20'
i.flow = i.flow[1:43]

# =========== a total of five kinds of weighting by topology ============= 
wt.0.nb =  c(rep(0,length(i.foc)),rep(1,length(i.flow)))
wt.0.f =  c(rep(1,length(i.foc)),rep(0,length(i.flow)))

wt.ab.nb = m.n.nb[c(i.foc,i.flow)]
wt.ab.f = c(m.n.foc[i.foc],m.n.nb[i.flow])

wt.a.nb = c(rep(0,length(i.foc)),m.n.nb[i.flow])
wt.b.nb = c(m.n.nb[i.foc],rep(0,length(i.flow)))
wt.b.f = c(m.n.foc[i.foc],rep(0,length(i.flow)))

# ===================== calculate prediction  ============================ 
# ** epx **
pred.ab.epx = apply(epx.surv,2,sum)/nrow(epx.surv)
pred.a.epx = apply(epx.surv*wt.0.nb,2,sum)/nrow(epx.surv)
pred.b.epx = apply(epx.surv*wt.0.f,2,sum)/nrow(epx.surv)

# ** mag **
pred.ab.nb.mag.x = abs(epx.surv)*wt.ab.nb
pred.ab.nb.mag = apply(pred.ab.nb.mag.x,2,sum)/nrow(epx.surv)
pred.ab.f.mag.x = abs(epx.surv)*wt.ab.f
pred.ab.f.mag = apply(pred.ab.f.mag.x,2,sum)/nrow(epx.surv)

pred.a.nb.mag.x = abs(epx.surv)*wt.a.nb
pred.a.nb.mag = apply(pred.a.nb.mag.x,2,sum)/nrow(epx.surv)
pred.b.nb.mag.x = abs(epx.surv)*wt.b.nb
pred.b.nb.mag = apply(pred.b.nb.mag.x,2,sum)/nrow(epx.surv)
pred.b.f.mag.x = abs(epx.surv)*wt.b.f
pred.b.f.mag = apply(pred.b.f.mag.x,2,sum)/nrow(epx.surv)

# ** spec **
pred.ab.nb.spec.1 = apply(abs(epx.surv),2,k.b.wt.dis.G,wt.ab.nb)/(nrow(epx.surv)-1)
pred.ab.nb.spec.2 = apply(abs(epx.surv),2,k.b.wt.dis.G.contra,wt.ab.nb)/(nrow(epx.surv)-1)
pred.ab.nb.spec.x = 0.5*(pred.ab.nb.spec.1 + pred.ab.nb.spec.2)
pred.ab.nb.spec = apply(pred.ab.nb.spec.x,2,sum)/nrow(epx.surv)

pred.ab.f.spec.1 = apply(abs(epx.surv),2,k.b.wt.dis.G,wt.ab.f)/(nrow(epx.surv)-1)
pred.ab.f.spec.2 = apply(abs(epx.surv),2,k.b.wt.dis.G.contra,wt.ab.f)/(nrow(epx.surv)-1)
pred.ab.f.spec.x = 0.5*(pred.ab.f.spec.1 + pred.ab.f.spec.2)
pred.ab.f.spec = apply(pred.ab.f.spec.x,2,sum)/nrow(epx.surv)

pred.a.nb.spec.1 = apply(abs(epx.surv),2,k.b.wt.dis.G,wt.a.nb)/(nrow(epx.surv)-1)
pred.a.nb.spec.2 = apply(abs(epx.surv),2,k.b.wt.dis.G.contra,wt.a.nb)/(nrow(epx.surv)-1)
pred.a.nb.spec.x = 0.5*(pred.a.nb.spec.1 + pred.a.nb.spec.2)
pred.a.nb.spec = apply(pred.a.nb.spec.x,2,sum)/nrow(epx.surv)

pred.b.nb.spec.1 = apply(abs(epx.surv),2,k.b.wt.dis.G,wt.b.nb)/(nrow(epx.surv)-1)
pred.b.nb.spec.2 = apply(abs(epx.surv),2,k.b.wt.dis.G.contra,wt.b.nb)/(nrow(epx.surv)-1)
pred.b.nb.spec.x = 0.5*(pred.b.nb.spec.1 + pred.b.nb.spec.2)
pred.b.nb.spec = apply(pred.b.nb.spec.x,2,sum)/nrow(epx.surv)

pred.b.f.spec.1 = apply(abs(epx.surv),2,k.b.wt.dis.G,wt.b.f)/(nrow(epx.surv)-1)
pred.b.f.spec.2 = apply(abs(epx.surv),2,k.b.wt.dis.G.contra,wt.b.f)/(nrow(epx.surv)-1)
pred.b.f.spec.x = 0.5*(pred.b.f.spec.1 + pred.b.f.spec.2)
pred.b.f.spec = apply(pred.b.f.spec.x,2,sum)/nrow(epx.surv)

# ** score **
pred.ab.nb.score = pred.ab.nb.mag/pred.ab.nb.spec
pred.ab.nb.score.g = pred.ab.nb.mag.x/pred.ab.nb.spec.x
pred.ab.f.score = pred.ab.f.mag/pred.ab.f.spec
pred.ab.f.score.g = pred.ab.f.mag.x/pred.ab.f.spec.x

pred.a.nb.score = pred.a.nb.mag/pred.a.nb.spec
pred.a.nb.score.g = pred.a.nb.mag.x/pred.a.nb.spec.x
pred.b.nb.score = pred.b.nb.mag/pred.b.nb.spec
pred.b.nb.score.g = pred.b.nb.mag.x/pred.b.nb.spec.x
pred.b.f.score = pred.b.f.mag/pred.b.f.spec
pred.b.f.score.g = pred.b.f.mag.x/pred.b.f.spec.x



# ================================================================================  	# ** ------ per gene plotting surv/hist/signal ---- **
# ------------------------------------------------------------------------------------------------------------------------------------------
library(geneplotter)
test.g = 'CBX5'
gind = match(test.g,rownames(epx.surv))

# **** =========== ab.nb =========== ****
gc.epx = epx.surv[gind,]
gc.mag = pred.ab.nb.mag.x[gind,]
gc.spec = pred.ab.nb.spec.x[gind,]
gc.score = pred.ab.nb.score.g[gind,]

# **** =========== ab.f =========== ****
gc.epx = epx.surv[gind,]
gc.mag = pred.ab.f.mag.x[gind,]
gc.spec = pred.ab.f.spec.x[gind,]
gc.score = pred.ab.f.score.g[gind,]

# **** =========== b.f =========== ****
gc.epx = epx.surv[gind,]
gc.mag = pred.b.f.mag.x[gind,]
gc.spec = pred.b.f.spec.x[gind,]
gc.score = pred.b.f.score.g[gind,]


# **** =========== Ratio (ab.nb) =========== ****
gc.epx = epx.surv[gind,]/pred.ab.epx
gc.mag = pred.ab.nb.mag.x[gind,]/pred.ab.nb.mag
gc.spec = pred.ab.nb.spec.x[gind,]/pred.ab.nb.spec
gc.score = pred.ab.nb.score.g[gind,]/pred.ab.nb.score


# **** =========== Ratio (ab.f) =========== ****
gc.epx = epx.surv[gind,]/pred.ab.epx
gc.mag = pred.ab.f.mag.x[gind,]/pred.ab.f.mag
gc.spec = pred.ab.f.spec.x[gind,]/pred.ab.f.spec
gc.score = pred.ab.f.score.g[gind,]/pred.ab.f.score


# **** =========== Ratio (b.f) =========== ****
gc.epx = epx.surv[gind,]/pred.b.epx
gc.mag = pred.b.f.mag.x[gind,]/pred.b.f.mag
gc.spec = pred.b.f.spec.x[gind,]/pred.b.f.spec
gc.score = pred.b.f.score.g[gind,]/pred.b.f.score


# **** =========== Ratio (a.nb) =========== ****
gc.epx = epx.surv[gind,]/pred.a.epx
gc.mag = pred.a.nb.mag.x[gind,]/pred.a.nb.mag
gc.spec = pred.a.nb.spec.x[gind,]/pred.a.nb.spec
gc.score = pred.a.nb.score.g[gind,]/pred.a.nb.score


# **** Global statistics **** # 
#w = cbind(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score)
#w = cbind(pred.ab.epx,pred.ab.f.mag,pred.ab.f.spec,pred.ab.f.score)
#w = cbind(pred.a.epx,pred.a.nb.mag,pred.a.nb.spec,pred.a.nb.score)
#w = cbind(pred.b.epx,pred.b.nb.mag,pred.b.nb.spec,pred.b.nb.score)
#w = cbind(pred.b.epx,pred.b.f.mag,pred.b.f.spec,pred.b.f.score)

# **** statistics **** # 
w = cbind(gc.epx,gc.mag,gc.spec,gc.score)
w0 = w[surv.st==0,]
w0.mns = surv.mns[surv.st==0]
w1 = w[surv.st==1,]
w1.mns = surv.mns[surv.st==1]
type = c('Exprs','Magnitude','Spectrum','SNR')


# * --------------------
# * Survival analysis
# * Histogram
# * --------------------
par(mfrow=c(2,2))
for (y in 1:4) {
#for (y in c(1,3)) {
#y = 3
xtitle = sprintf("%s - %s",test.g,type[y])
#xtitle = sprintf("Hubs-%s",type[y])
#plot.surv(w[,y],surv.mns,surv.st,xtitle)
#}

br = seq(floor(min(w[,y])),ceiling(max(w[,y])),length=35)
br = round(br,digit=1)
cols <- c('lightgrey','pink')
colss <- c('pink','lightgrey')
x0 = w0[,y]
x1 = w1[,y]
xpv = wilcox.test(x0,x1)$p.value
txt= sprintf("Wilcox.test Pv=%.3f", round(xpv,digit=3))
histStack(list(x1,x0), breaks=br, col=colss,cex.axis=0.9)
#legend(15,80,c(xtitle,"alive","dead",txt),pch=c(-1,7,7,-1),col=c("",cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")
legend(10,55,c(xtitle,"alive","dead",txt),pch=c(-1,7,7,-1),col=c("",cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")
}


# * ----------
# * dynamics
# * ----------
par(mfrow=c(2,2))
#y = 1
for (y in 1:4) {
#for (y in c(1,3)) {
xtitle = sprintf("%s - %s",test.g,type[y])
ind.x0 = order(w0.mns)
w0.x0=w0[ind.x0,]
x0.axis = w0.mns[ind.x0]
ind.x1 = order(w1.mns)
w1.x1=w1[ind.x1,]
x1.axis = w1.mns[ind.x1]

plot(w0.x0[,y]~x0.axis,pch=20,cex=0.4,col='black',xlim=c(0,max(surv.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='OS (months)',ylab=xtitle)
lines(w0.x0[,y]~x0.axis,col='lightgray',xlim=c(0,max(surv.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='OS (months)',ylab=xtitle)
points(w1.x1[,y]~x1.axis,pch=20,cex=0.4,col='red',xlim=c(0,max(surv.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='OS (months)',ylab=xtitle)
lines(w1.x1[,y]~x1.axis,col='pink',xlim=c(0,max(surv.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='OS (months)',ylab=xtitle)

ks.ck.pv = unlist(ks.boot(w0.x0[,y],w1.x1[,y])$ks[2])
pvv = sprintf("P-value (K-S) = %.3f", round(ks.ck.pv,digit=3))
legend(80,max(w[,y])+1,c("alive","dead",pvv),pch=c(20,20,-1),lty=c(1,1,-1),col=c(cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")
}



# =================================================================================================================================================
# corr surv time
# =================================================================================================================================================

spcol.0 <- colorRampPalette(c("yellow", "orange"))(nrow(w0))						
spcol.1=rep('red',nrow(w1))															

ind.ex = rep(T,length(surv.mns))					# * ---- remove surv = 0
ind.ex[surv.mns==0] = F
gc.mns = surv.mns[ind.ex]
gc.st = surv.st[ind.ex]
w = w[ind.ex,]
w0 = w[gc.st==0,]
w0.mns = gc.mns[gc.st==0]
w1 = w[gc.st==1,]
w1.mns = gc.mns[gc.st==1]
spcol.0 <- colorRampPalette(c("yellow", "orange"))(nrow(w0))						
spcol.1=rep('red',nrow(w1))															

# ****
#k = y
# ****
par(mfrow=c(2,2))
for (k in 1:4) {
#for (y in c(1,3)) {

plot(w0.mns~w0[,k],col=spcol.0,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(gc.mns)),ylab="OS(months)",xlab= type[k])
abline(glm(w0.mns~w0[,k],family=gaussian),lty=2,lwd=1.8,col='orange')
stat = summary(glm(w0.mns~w0[,k],family=gaussian))
e.pv.0 = sprintf("Alive: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
c.stat = cor.test(w0.mns,w0[,k])
e.corr.0 = sprintf("Corr=%.2f Pv=%.2f",c.stat[4],c.stat[3])
#legend(median(w[,k]),max(gc.mns),c("meta-free","meta",e.pv.0),cex=0.6,col=c('orange','red'),pch=c(1,1,-1),text.col ='black', bg ='white',bty="n")
legend(min(w[,k]),max(gc.mns),c(e.pv.0),cex=0.6,text.col ='black', bg ='white',bty="n")
cat(e.corr.0,"\n")

points(w1.mns~w1[,k],col=spcol.1,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(surv.mns)),ylab="OS (months)",xlab= type[k])
abline(glm(w1.mns~w1[,k],family=gaussian),lty=2,lwd=1.8,col='red')
stat = summary(glm(w1.mns~w1[,k],family=gaussian))
e.pv.1 = sprintf("Dead: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
c.stat.1 = cor.test(w1.mns,w1[,k])
e.corr.1 = sprintf("Corr=%.2f Pv=%.2f",c.stat.1[4],c.stat.1[3])
legend(min(w[,k]),max(gc.mns)-5,c(e.pv.1),cex=0.6,text.col ='black', bg ='white',bty="n")
cat(e.corr.1,"\n")

abline(glm(gc.mns~w[,k],family=gaussian),lty=2,lwd=1.8,col='black')
stat = summary(glm(gc.mns~w[,k],family=gaussian))
e.pv.01 = sprintf("All: E=%.2f Pv= %.2f",stat[[12]][2],stat[[12]][8])
c.stat.01 = cor.test(gc.mns,w[,k])
e.corr.01 = sprintf("Corr=%.2f Pv=%.2f",c.stat.01[4],c.stat.01[3])
legend(min(w[,k]),max(gc.mns)-10,c(e.pv.01),cex=0.6,text.col ='black', bg ='white',bty="n")
cat(e.corr.01,"\n\n\n")

}

	

	
	
# =================================================================================================================================================
# corr 1/surv(T)
# =================================================================================================================================================	
# 
#k=1
par(mfrow=c(2,2))
for (k in 1:4) {

plot(1/w0.mns~w0[,k],col=spcol.0,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(1/gc.mns)),ylab="1 / MFS(months)",xlab= type[k])
abline(glm(1/w0.mns~w0[,k],family=gaussian),lty=2,lwd=1.8,col='orange')
stat = summary(glm(1/w0.mns~w0[,k],family=gaussian))
e.pv.0 = sprintf("Alive: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
c.stat = cor.test(1/w0.mns,w0[,k])
e.corr.0 = sprintf("Corr=%.2f Pv=%.2f",c.stat[4],c.stat[3])
legend(min(w[,k]),max(1/gc.mns),c("meta-free","meta",e.pv.0),cex=0.6,col=c('orange','red'),pch=c(1,1,-1),text.col ='black', bg ='white',bty="n")
cat("Alive",e.corr.0,"\n")
	
points(1/w1.mns~w1[,k],col=spcol.1,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(surv.mns)),ylab="MFS (months)",xlab= type[k])
abline(glm(1/w1.mns~w1[,k],family=gaussian),lty=2,lwd=1.8,col='red')
stat = summary(glm(1/w1.mns~w1[,k],family=gaussian))
e.pv.1 = sprintf("Dead: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
c.stat.1 = cor.test(1/w1.mns,w1[,k])
e.corr.1 = sprintf("Corr=%.2f Pv=%.2f",c.stat.1[4],c.stat.1[3])
legend(min(w[,k]),max(1/gc.mns)-0.1,c(e.pv.1),cex=0.6,text.col ='black', bg ='white',bty="n")
cat("Dead",e.corr.1,"\n")

abline(glm(1/gc.mns~w[,k],family=gaussian),lty=2,lwd=1.8,col='black')
stat = summary(glm(1/gc.mns~w[,k],family=gaussian))
e.pv.01 = sprintf("All: E=%.2f Pv= %.2f",stat[[12]][2],stat[[12]][8])
c.stat.01 = cor.test(gc.mns,w[,k])
e.corr.01 = sprintf("Corr=%.2f Pv=%.2f",c.stat.01[4],c.stat.01[3])
legend(min(w[,k]),max(1/gc.mns)-0.15,c(e.pv.01),cex=0.6,text.col ='black', bg ='white',bty="n")
cat("All",e.corr.01,"\n\n")

}






