# ----------------------------------------------------------------------------------------------------------------------------------------- #
# * meta 
# * merge all datasets
# ------------------------------------------------------------------------------------------------------------------------------------------
library('survival')
library('affyQCReport')
library('limma')
library('Mfuzz')
library('fdrtool')
library('Matching')
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

file = c('msk1','msk2','smc8894')

dt.files = vector('list',3)

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

# --------------------------------------------------------------- # 
n=3
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
# *** meta *** #
epx.met = e.sym
met.st = c()
met.mns = c()
for (e in 1:2) {
	dd = dt.files[[e]]$data
	dd.1 = dd[,match(e.sym,colnames(dd),nomatch=0)]
	dd.2 = cbind(dd.1,rep(0,nrow(dd.1)))
	colnames(dd.2) = e.sym
	epx.met = rbind(epx.met,dd.2)
	met.st = append(met.st,dt.files[[e]]$met.st)
	met.mns = append(met.mns,dt.files[[e]]$met.mns)
}
e = 3 
dd = dt.files[[e]]$data
dd.1 = dd[,match(e.sym,colnames(dd),nomatch=0)]
epx.met = rbind(epx.met,dd.1)
met.st = append(met.st,dt.files[[e]]$surv.st)
met.mns = append(met.mns,dt.files[[e]]$surv.mns)
epx.met = epx.met[2:nrow(epx.met),]
epx.met = t(epx.met)
k = dim(epx.met)
epx.met = as.numeric(epx.met)
dim(epx.met) = k
rownames(epx.met) = e.sym

	
# ******************************************************************************************************************************* #
epx.met = epx.met[1:68,]										# remove missing 'DDX20' ;; USE ONLY HUB GENES !! #
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
pred.ab.epx = apply(epx.met,2,sum)/nrow(epx.met)
pred.a.epx = apply(epx.met*wt.0.nb,2,sum)/nrow(epx.met)
pred.b.epx = apply(epx.met*wt.0.f,2,sum)/nrow(epx.met)

# ** mag **
pred.ab.nb.mag.x = abs(epx.met)*wt.ab.nb
pred.ab.nb.mag = apply(pred.ab.nb.mag.x,2,sum)/nrow(epx.met)
pred.ab.f.mag.x = abs(epx.met)*wt.ab.f
pred.ab.f.mag = apply(pred.ab.f.mag.x,2,sum)/nrow(epx.met)

pred.a.nb.mag.x = abs(epx.met)*wt.a.nb
pred.a.nb.mag = apply(pred.a.nb.mag.x,2,sum)/nrow(epx.met)
pred.b.nb.mag.x = abs(epx.met)*wt.b.nb
pred.b.nb.mag = apply(pred.b.nb.mag.x,2,sum)/nrow(epx.met)
pred.b.f.mag.x = abs(epx.met)*wt.b.f
pred.b.f.mag = apply(pred.b.f.mag.x,2,sum)/nrow(epx.met)

# ** spec **
pred.ab.nb.spec.1 = apply(abs(epx.met),2,k.b.wt.dis.G,wt.ab.nb)/(nrow(epx.met)-1)
pred.ab.nb.spec.2 = apply(abs(epx.met),2,k.b.wt.dis.G.contra,wt.ab.nb)/(nrow(epx.met)-1)
pred.ab.nb.spec.x = 0.5*(pred.ab.nb.spec.1 + pred.ab.nb.spec.2)
pred.ab.nb.spec = apply(pred.ab.nb.spec.x,2,sum)/nrow(epx.met)

pred.ab.f.spec.1 = apply(abs(epx.met),2,k.b.wt.dis.G,wt.ab.f)/(nrow(epx.met)-1)
pred.ab.f.spec.2 = apply(abs(epx.met),2,k.b.wt.dis.G.contra,wt.ab.f)/(nrow(epx.met)-1)
pred.ab.f.spec.x = 0.5*(pred.ab.f.spec.1 + pred.ab.f.spec.2)
pred.ab.f.spec = apply(pred.ab.f.spec.x,2,sum)/nrow(epx.met)

pred.a.nb.spec.1 = apply(abs(epx.met),2,k.b.wt.dis.G,wt.a.nb)/(nrow(epx.met)-1)
pred.a.nb.spec.2 = apply(abs(epx.met),2,k.b.wt.dis.G.contra,wt.a.nb)/(nrow(epx.met)-1)
pred.a.nb.spec.x = 0.5*(pred.a.nb.spec.1 + pred.a.nb.spec.2)
pred.a.nb.spec = apply(pred.a.nb.spec.x,2,sum)/nrow(epx.met)

pred.b.nb.spec.1 = apply(abs(epx.met),2,k.b.wt.dis.G,wt.b.nb)/(nrow(epx.met)-1)
pred.b.nb.spec.2 = apply(abs(epx.met),2,k.b.wt.dis.G.contra,wt.b.nb)/(nrow(epx.met)-1)
pred.b.nb.spec.x = 0.5*(pred.b.nb.spec.1 + pred.b.nb.spec.2)
pred.b.nb.spec = apply(pred.b.nb.spec.x,2,sum)/nrow(epx.met)

pred.b.f.spec.1 = apply(abs(epx.met),2,k.b.wt.dis.G,wt.b.f)/(nrow(epx.met)-1)
pred.b.f.spec.2 = apply(abs(epx.met),2,k.b.wt.dis.G.contra,wt.b.f)/(nrow(epx.met)-1)
pred.b.f.spec.x = 0.5*(pred.b.f.spec.1 + pred.b.f.spec.2)
pred.b.f.spec = apply(pred.b.f.spec.x,2,sum)/nrow(epx.met)

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


# ================================================================================  			# ** ------ surv T or 1/T GLM model --- ** #

# Continue to  proc_meta_g_Corr.r

# ------------------------------------------------------------------------------------------------------------------------------------------
par(mfrow=c(2,2))

k=4
bin = seq(6,140,by=6)
type = c('Ave. Exprs','Magnitude','Spectrum','SNR')
# ======================= #
w = cbind(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score)
w0 = w[met.st==0,]
w0.mns = met.mns[met.st==0]
w1 = w[met.st==1,]
w1.mns = met.mns[met.st==1]
spcol.0 <- colorRampPalette(c("yellow", "orange"))(nrow(w0))						
spcol.1=rep('red',nrow(w1))															

#for (k in 1:2) {
	plot(w0.mns~w0[,k],col=spcol.0,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(met.mns)),ylab="MFS (months)",xlab= type[k])
	abline(glm(w0.mns~w0[,k],family=gaussian),lty=2,lwd=1.8,col='darkgrey')
	stat = summary(glm(w0.mns~w0[,k],family=gaussian))
	e.pv.0 = sprintf("Meta-free: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(met.mns),c("meta-free","meta",e.pv.0),cex=0.6,col=c('orange','red'),pch=c(1,1,-1),text.col ='black', bg ='white',bty="n")
	
	points(w1.mns~w1[,k],col=spcol.1,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(met.mns)),ylab="MFS (months)",xlab= type[k])
	abline(glm(w1.mns~w1[,k],family=gaussian),lty=2,lwd=1.8,col='brown')
	stat = summary(glm(w1.mns~w1[,k],family=gaussian))
	e.pv.1 = sprintf("Meta: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(met.mns)-10,c(e.pv.1),cex=0.6,text.col ='black', bg ='white',bty="n")

	abline(glm(met.mns~w[,k],family=gaussian),lty=2,lwd=1.8,col='blue')
	stat = summary(glm(met.mns~w[,k],family=gaussian))
	e.pv.01 = sprintf("All: E=%.2f Pv= %.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(met.mns)-13,c(e.pv.01),cex=0.6,text.col ='black', bg ='white',bty="n")
#}

# ============================
ind.ex = rep(T,length(met.mns))					# * ---- remove surv = 0
ind.ex[met.mns==0] = F
gc.mns = met.mns[ind.ex]
gc.st = met.st[ind.ex]
w = w[ind.ex,]
w0 = w[gc.st==0,]
w0.mns = gc.mns[gc.st==0]
w1 = w[gc.st==1,]
w1.mns = gc.mns[gc.st==1]
spcol.0 <- colorRampPalette(c("yellow", "orange"))(nrow(w0))						
spcol.1=rep('red',nrow(w1))															

#for (k in 3:4) {
	plot(1/w0.mns~w0[,k],col=spcol.0,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(1/gc.mns)),ylab="1 / MFS(months)",xlab= type[k])
	abline(glm(1/w0.mns~w0[,k],family=gaussian),lty=2,lwd=1.8,col='darkgrey')
	stat = summary(glm(1/w0.mns~w0[,k],family=gaussian))
	e.pv.0 = sprintf("Meta-free: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(1/gc.mns),c("meta-free","meta",e.pv.0),cex=0.6,col=c('orange','red'),pch=c(1,1,-1),text.col ='black', bg ='white',bty="n")
	
	points(1/w1.mns~w1[,k],col=spcol.1,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(met.mns)),ylab="MFS (months)",xlab= type[k])
	abline(glm(1/w1.mns~w1[,k],family=gaussian),lty=2,lwd=1.8,col='brown')
	stat = summary(glm(1/w1.mns~w1[,k],family=gaussian))
	e.pv.1 = sprintf("Meta: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(1/gc.mns)-0.08,c(e.pv.1),cex=0.6,text.col ='black', bg ='white',bty="n")

	abline(glm(1/gc.mns~w[,k],family=gaussian),lty=2,lwd=1.8,col='blue')
	stat = summary(glm(1/gc.mns~w[,k],family=gaussian))
	e.pv.01 = sprintf("All: E=%.2f Pv= %.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(1/gc.mns)-0.1,c(e.pv.01),cex=0.6,text.col ='black', bg ='white',bty="n")
#}




# ================================================================================  		# ** ------ Wilcox.test (box/histStack) --- ** #
# ------------------------------------------------------------------------------------------------------------------------------------------
library(geneplotter)

type = c('Ave. Exprs','Magnitude','Spectrum','SNR')
w = cbind(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score)
w.ind = met.st
w0 = w[met.st==0,]
w0.mns = met.mns[met.st==0]
w1 = w[met.st==1,]
w1.mns = met.mns[met.st==1]

# ======================= #

k=4
xtitle = type[k]
#for (k in 1:4) {
	br = seq(floor(min(w[,k])),ceiling(max(w[,k])),length=30)
	br = round(br,digit=1)
	cols <- c('lightgrey','pink')
	colss <- c('pink','lightgrey')
	x0 = w0[,k]
	x1 = w1[,k]
	
	xpv = wilcox.test(x0,x1)$p.value
	txt= sprintf("Wilcox.test Pv=%.2f", xpv)
	histStack(list(x1,x0), breaks=br, col=colss,cex.axis=0.9)
	legend(20,65,c(xtitle,"meta-free","meta",txt),pch=c(-1,7,7,-1),col=c("",cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")
#}

# ======================= dynamics --------------- #
ind.x0 = order(w0.mns)
w0.x0=w0[ind.x0,]
x0.axis = w0.mns[ind.x0]
ind.x1 = order(w1.mns)
w1.x1=w1[ind.x1,]
x1.axis = w1.mns[ind.x1]

plot(w0.x0[,k]~x0.axis,pch=20,cex=0.4,col='black',xlim=c(0,max(met.mns)),ylim=c((min(w[,k])-1),max(w[,k])+1),xlab='MFS (months)',ylab=xtitle)
lines(w0.x0[,k]~x0.axis,col='lightgray',xlim=c(0,max(met.mns)),ylim=c((min(w[,k])-1),max(w[,k])+1),xlab='MFS (months)',ylab=xtitle)
points(w1.x1[,k]~x1.axis,pch=20,cex=0.4,col='red',xlim=c(0,max(met.mns)),ylim=c((min(w[,k])-1),max(w[,k])+1),xlab='MFS (months)',ylab=xtitle)
lines(w1.x1[,k]~x1.axis,col='pink',xlim=c(0,max(met.mns)),ylim=c((min(w[,k])-1),max(w[,k])+1),xlab='MFS (months)',ylab=xtitle)

ks.ck.pv = unlist(ks.boot(w0.x0[,k],w1.x1[,k])$ks[2])
pvv = sprintf("P-value (K-S) = %.3f", round(ks.ck.pv,digit=3))
legend(80,max(w[,k])+1,c("meta-free","meta",pvv),pch=c(20,20,-1),lty=c(1,1,-1),col=c(cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")
#}




# ================================================================================  	# ** ------ per gene plotting surv/hist/signal ---- **
# ------------------------------------------------------------------------------------------------------------------------------------------
library(geneplotter)
test.g = 'CBX5'
gind = match(test.g,rownames(epx.met))

gc.epx = epx.met[gind,]
gc.mag = pred.ab.nb.mag.x[gind,]
gc.spec = pred.ab.nb.spec.x[gind,]
gc.score = pred.ab.nb.score.g[gind,]
w = cbind(gc.epx,gc.mag,gc.spec,gc.score)
w0 = w[met.st==0,]
w0.mns = met.mns[met.st==0]
w1 = w[met.st==1,]
w1.mns = met.mns[met.st==1]

par(mfrow=c(3,3))
type = c('Exprs','Magnitude','Spectrum','SNR')
#for (y in c(1,3,4)) {
y = 1
xtitle = sprintf("%s - %s",test.g,type[y])
plot.surv(w[,y],met.mns,met.st,xtitle)

br = seq(floor(min(w[,y])),ceiling(max(w[,y])),length=30)
br = round(br,digit=1)
cols <- c('lightgrey','pink')
colss <- c('pink','lightgrey')
x0 = w0[,y]
x1 = w1[,y]
xpv = wilcox.test(x0,x1)$p.value
txt= sprintf("Wilcox.test Pv=%.3f", round(xpv,digit=3))
histStack(list(x1,x0), breaks=br, col=colss,cex.axis=0.9)
#legend(15,80,c(xtitle,"meta-free","meta",txt),pch=c(-1,7,7,-1),col=c("",cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")
legend(10,55,c(xtitle,"meta-free","meta",txt),pch=c(-1,7,7,-1),col=c("",cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")

# * ----------
# * dynamics
# * ----------
par(mfrow=c(3,1))
#y = 1
for (y in c(1,3,4)) {
xtitle = sprintf("%s - %s",test.g,type[y])
ind.x0 = order(w0.mns)
w0.x0=w0[ind.x0,]
x0.axis = w0.mns[ind.x0]
ind.x1 = order(w1.mns)
w1.x1=w1[ind.x1,]
x1.axis = w1.mns[ind.x1]

plot(w0.x0[,y]~x0.axis,pch=20,cex=0.4,col='black',xlim=c(0,max(met.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='MFS (months)',ylab=xtitle)
lines(w0.x0[,y]~x0.axis,col='lightgray',xlim=c(0,max(met.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='MFS (months)',ylab=xtitle)
points(w1.x1[,y]~x1.axis,pch=20,cex=0.4,col='red',xlim=c(0,max(met.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='MFS (months)',ylab=xtitle)
lines(w1.x1[,y]~x1.axis,col='pink',xlim=c(0,max(met.mns)),ylim=c((min(w[,y])-1),max(w[,y])+1),xlab='MFS (months)',ylab=xtitle)

ks.ck.pv = unlist(ks.boot(w0.x0[,y],w1.x1[,y])$ks[2])
pvv = sprintf("P-value (K-S) = %.3f", round(ks.ck.pv,digit=3))
legend(80,max(w[,y])+1,c("meta-free","meta",pvv),pch=c(20,20,-1),lty=c(1,1,-1),col=c(cols,""),cex=0.7,text.col = 'black', bg ='white',bty="n")
}


# **************************
# corr with 1/Surv(T)
# ==========================

k=3
test.g = 'DNMT1'
gind = match(test.g,rownames(epx.met))

gc.epx = epx.met[gind,]
gc.mag = pred.ab.nb.mag.x[gind,]
gc.spec = pred.ab.nb.spec.x[gind,]
gc.score = pred.ab.nb.score.g[gind,]
w = cbind(gc.epx,gc.mag,gc.spec,gc.score)

type = c('Ave. Exprs','Magnitude','Spectrum','SNR')
spcol.0 <- colorRampPalette(c("yellow", "orange"))(nrow(w0))						
spcol.1=rep('red',nrow(w1))															

ind.ex = rep(T,length(met.mns))					# * ---- remove surv = 0
ind.ex[met.mns==0] = F
gc.mns = met.mns[ind.ex]
gc.st = met.st[ind.ex]
w = w[ind.ex,]
w0 = w[gc.st==0,]
w0.mns = gc.mns[gc.st==0]
w1 = w[gc.st==1,]
w1.mns = gc.mns[gc.st==1]
spcol.0 <- colorRampPalette(c("yellow", "orange"))(nrow(w0))						
spcol.1=rep('red',nrow(w1))															

	plot(1/w0.mns~w0[,k],col=spcol.0,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(1/gc.mns)),ylab="1 / MFS(months)",xlab= type[k])
	abline(glm(1/w0.mns~w0[,k],family=gaussian),lty=2,lwd=1.8,col='darkgrey')
	stat = summary(glm(1/w0.mns~w0[,k],family=gaussian))
	e.pv.0 = sprintf("Meta-free: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(1/gc.mns),c("meta-free","meta",e.pv.0),cex=0.6,col=c('orange','red'),pch=c(1,1,-1),text.col ='black', bg ='white',bty="n")
	
	points(1/w1.mns~w1[,k],col=spcol.1,xlim=c(min(w[,k]),max(w[,k])),ylim=c(0,max(met.mns)),ylab="MFS (months)",xlab= type[k])
	abline(glm(1/w1.mns~w1[,k],family=gaussian),lty=2,lwd=1.8,col='brown')
	stat = summary(glm(1/w1.mns~w1[,k],family=gaussian))
	e.pv.1 = sprintf("Meta: E=%.2f Pv=%.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(1/gc.mns)-0.08,c(e.pv.1),cex=0.6,text.col ='black', bg ='white',bty="n")

	abline(glm(1/gc.mns~w[,k],family=gaussian),lty=2,lwd=1.8,col='blue')
	stat = summary(glm(1/gc.mns~w[,k],family=gaussian))
	e.pv.01 = sprintf("All: E=%.2f Pv= %.2f",stat[[12]][2],stat[[12]][8])
	legend(min(w[,k]),max(1/gc.mns)-0.1,c(e.pv.01),cex=0.6,text.col ='black', bg ='white',bty="n")

