# ----------------------------------------------------------------------------------------------------------------------------------------- #
# * meta
# * lvcom 
# ------------------------------------------------------------------------------------------------------------------------------------------
library('survival')
library('affyQCReport')
library('limma')
library('Mfuzz')
#library('fdrtool')
source("sub_func.r")
source("method.r")
source("method_surv.r")
source('proc_g_plot.r')


# -------- probe-Lung-specific ------ # 																					# --- annotation --- #
mark <-  read.delim("lvcom_cnvt.txt", sep="\t",header=T)
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
	
	inf = sprintf("lvcom_nwEpx_%s.txt",file[n])
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
	
	inf = sprintf("lvcom_nwEpx_%s.txt",file[n])
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
	#dd.2 = cbind(dd.1,rep(0,nrow(dd.1)))
	dd.2=dd.1
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

# ** wt.exprs **
pred.ab.nb.wtEpx.x = epx.met*wt.ab.nb
pred.ab.nb.wtEpx = apply(pred.ab.nb.wtEpx.x,2,sum)/nrow(epx.met)
pred.ab.f.wtEpx.x = epx.met*wt.ab.f
pred.ab.f.wtEpx = apply(pred.ab.f.wtEpx.x,2,sum)/nrow(epx.met)

pred.a.nb.wtEpx.x = epx.met*wt.a.nb
pred.a.nb.wtEpx = apply(pred.a.nb.wtEpx.x,2,sum)/nrow(epx.met)
pred.b.nb.wtEpx.x = epx.met*wt.b.nb
pred.b.nb.wtEpx = apply(pred.b.nb.wtEpx.x,2,sum)/nrow(epx.met)
pred.b.f.wtEpx.x = epx.met*wt.b.f
pred.b.f.wtEpx = apply(pred.b.f.wtEpx.x,2,sum)/nrow(epx.met)

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


# ================================================================================  				# ** ------ 1 test survival --- ** #
test.ab.nb = report.nw.surv(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score,met.mns,met.st) 
test.ab.f = report.nw.surv(pred.ab.epx,pred.ab.f.mag,pred.ab.f.spec,pred.ab.f.score,met.mns,met.st) 
test.a.nb = report.nw.surv(pred.a.epx,pred.a.nb.mag,pred.a.nb.spec,pred.a.nb.score,met.mns,met.st) 
test.b.nb = report.nw.surv(pred.b.epx,pred.b.nb.mag,pred.b.nb.spec,pred.b.nb.score,met.mns,met.st) 
test.b.f = report.nw.surv(pred.b.epx,pred.b.f.mag,pred.b.f.spec,pred.b.f.score,met.mns,met.st) 

test = c('epx','mag','spec','snr')
test = data.matrix(cbind(test,test.ab.nb,test.ab.f,test.a.nb,test.b.nb,test.b.f))
write.table(test,"lvcom_meta_nw.txt",sep="\t",col.names =T, row.names=F, quote=F)

# == use wt.Exprs ==
test.ab.nb.2 = report.nw.surv(pred.ab.epx,pred.ab.nb.wtEpx,pred.ab.nb.spec,pred.ab.nb.score,met.mns,met.st)
test.ab.f.2 = report.nw.surv(pred.ab.epx,pred.ab.f.wtEpx,pred.ab.f.spec,pred.ab.f.score,met.mns,met.st) 
test.a.nb.2 = report.nw.surv(pred.a.epx,pred.a.nb.wtEpx,pred.a.nb.spec,pred.a.nb.score,met.mns,met.st)
test.b.nb.2 = report.nw.surv(pred.b.epx,pred.b.nb.wtEpx,pred.b.nb.spec,pred.b.nb.score,met.mns,met.st)
test.b.f.2 = report.nw.surv(pred.b.epx,pred.b.f.wtEpx,pred.b.f.spec,pred.b.f.score,met.mns,met.st)

test.2 = c('epx','wtEpx','spec','snr')
test = data.matrix(cbind(test.2,test.ab.nb.2,test.ab.f.2,test.a.nb.2,test.b.nb.2,test.b.f.2))
write.table(test,"lvcom_meta_nw_2.txt",sep="\t",col.names =T, row.names=F, quote=F)


# == *** plotting *** ==																			*** ------ all genes (nw) plotting --- **
plot.surv(pred.a.nb.score,met.mns,met.st,'merged metastasis-free survivals by averaged exprs')

plot.dyn.nw(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score,met.mns,met.st) 


# == * corr - 1/Surv(T) * ==	------------------------------------------------------------------	* corr - 1/Surv(T) * =============
ind.ex = rep(T,length(met.mns))					# * ---- remove surv = 0
ind.ex[met.mns==0] = F
gc.mns = met.mns[ind.ex]
gc.st = met.st[ind.ex]
gc.epx = pred.ab.epx[ind.ex]
gc.mag = pred.ab.nb.mag[ind.ex]
gc.spec = pred.ab.nb.spec[ind.ex]
gc.score = pred.ab.nb.score[ind.ex]

ind.q = order(gc.epx)
ind.st = gc.st[ind.q]
spcol.2 <- colorRampPalette(c("yellow", "orange"))(length(gc.spec))						
spcol.2[ind.st==1] = 'red'															

par(mfrow=c(3,5))

# ------------------------------------------------------------------------------------------------- * averaged expression * ------- **
plot(1/(gc.mns[ind.q])~gc.epx[ind.q],col=spcol.2,ylab="1 / MFS (months)",xlab= "Averaged expression")
abline(glm(1/(gc.mns)~gc.epx,family=gaussian),lty=2,lwd=1.5,col='black')
stat = summary(glm(1/(gc.mns)~gc.epx,family=gaussian),lty=2,lwd=1.5,col='black')
txt.esti = sprintf("Estimate (GLM) = %.3f",stat[[12]][2])
txt.pv = sprintf("P-value = %.3f",stat[[12]][8])
legend(min(gc.epx),max(1/gc.mns),c("meta-free","meta",txt.esti,txt.pv),cex=0.6,col=c('orange','red'),pch=c(1,1,-1,-1),text.col ='black', bg ='white',bty="n")

MFS = gc.mns
x = gc.epx
p.sort = summary(x)														  
Q = rep(0,length(x))
Q[x<=p.sort[[2]]] = 1
Q[x>p.sort[[2]]] = 2
Q[x>p.sort[[3]]] = 3
Q[x>p.sort[[5]]] = 4
hist(MFS[Q==1],col='darkgrey',breaks=30,main = 'Exprs - Q1')
hist(MFS[Q==2],col=2,breaks=30,main = 'Exprs - Q2')
hist(MFS[Q==3],col=3,breaks=30,main = 'Exprs - Q3')
hist(MFS[Q==4],col=4,breaks=30,main = 'Exprs - Q4')
# ---------------------------------------------------------------------------------------------------- * spec * ------------------ **
plot(1/(gc.mns[ind.q])~gc.spec[ind.q],col=spcol.2,ylab="1 / MFS (months)",xlab= "Spectrum")
abline(glm(1/(gc.mns)~gc.spec,family=gaussian),lty=2,lwd=1.5,col='black')
stat = summary(glm(1/(gc.mns)~gc.spec,family=gaussian),lty=2,lwd=1.5,col='black')
txt.esti = sprintf("Estimate (GLM) = %.3f",stat[[12]][2])
txt.pv = sprintf("P-value = %.3f",stat[[12]][8])
legend(min(gc.spec),max(1/gc.mns),c("meta-free","meta",txt.esti,txt.pv),cex=0.6,col=c('orange','red'),pch=c(1,1,-1,-1),text.col ='black', bg ='white',bty="n")

x = gc.spec
p.sort = summary(x)														  
Q = rep(0,length(x))
Q[x<=p.sort[[2]]] = 1
Q[x>p.sort[[2]]] = 2
Q[x>p.sort[[3]]] = 3
Q[x>p.sort[[5]]] = 4
hist(MFS[Q==1],col='darkgrey',breaks=30,main = 'Spectrum - Q1')
hist(MFS[Q==2],col=2,breaks=30,main = 'Spectrum - Q2')
hist(MFS[Q==3],col=3,breaks=30,main = 'Spectrum - Q3')
hist(MFS[Q==4],col=4,breaks=30,main = 'Spectrum - Q4')
# ---------------------------------------------------------------------------------------------------- * SNR * ------------------ **
plot(1/(gc.mns[ind.q])~gc.score[ind.q],col=spcol.2,ylab="1 / MFS (months)",xlab= "SNR")
abline(glm(1/(gc.mns)~gc.score,family=gaussian),lty=2,lwd=1.5,col='black')
stat = summary(glm(1/(gc.mns)~gc.score,family=gaussian),lty=2,lwd=1.5,col='black')
txt.esti = sprintf("Estimate (GLM) = %.3f",stat[[12]][2])
txt.pv = sprintf("P-value = %.3f",stat[[12]][8])
legend(min(gc.score),max(1/gc.mns),c("meta-free","meta",txt.esti,txt.pv),cex=0.6,col=c('orange','red'),pch=c(1,1,-1,-1),text.col ='black', bg ='white',bty="n")

x = gc.score
p.sort = summary(x)														  
Q = rep(0,length(x))
Q[x<=p.sort[[2]]] = 1
Q[x>p.sort[[2]]] = 2
Q[x>p.sort[[3]]] = 3
Q[x>p.sort[[5]]] = 4
hist(MFS[Q==1],col='darkgrey',breaks=30,main = 'SNR - Q1')
hist(MFS[Q==2],col=2,breaks=30,main = 'SNR - Q2')
hist(MFS[Q==3],col=3,breaks=30,main = 'SNR - Q3')
hist(MFS[Q==4],col=4,breaks=30,main = 'SNR - Q4')


# ================================================================================  	# ** ------ 2 per gene all data test survival --- ** =====
test.g = rownames(epx.met)
tg.ab.nb = report.g.surv(test.g,epx.met,pred.ab.nb.mag.x,pred.ab.nb.spec.x,pred.ab.nb.score.g,met.mns,met.st,"z_all_meta_g_abnb.txt") 
tg.ab.f = report.g.surv(test.g,epx.met,pred.ab.f.mag.x,pred.ab.f.spec.x,pred.ab.f.score.g,met.mns,met.st,"z_all_meta_g_abf.txt") 
test.g = m.sym[i.flow]
tg.a.nb = report.g.surv(test.g,epx.met,pred.a.nb.mag.x,pred.a.nb.spec.x,pred.a.nb.score.g,met.mns,met.st,"z_all_meta_g_a_nb.txt") 
test.g = m.sym[i.foc]
tg.b.nb = report.g.surv(test.g,epx.met,pred.b.nb.mag.x,pred.b.nb.spec.x,pred.b.nb.score.g,met.mns,met.st,"z_all_meta_g_b_nb.txt") 
tg.b.f = report.g.surv(test.g,epx.met,pred.b.f.mag.x,pred.b.f.spec.x,pred.b.f.score.g,met.mns,met.st,"z_all_meta_g_b_f.txt") 


# ================================================================================  	# ** ----- 3 per gene all data cor.test 1 / T(surv) --- **
test.g = rownames(epx.met)
tg.ab.nb = report.g.cor.invT(test.g,epx.met,pred.ab.nb.mag.x,pred.ab.nb.spec.x,pred.ab.nb.score.g,met.mns,"z_allmeta_corinvT_g_abnb.txt") 
tg.ab.f = report.g.cor.invT(test.g,epx.met,pred.ab.f.mag.x,pred.ab.f.spec.x,pred.ab.f.score.g,met.mns,"z_allmeta_corinvT_g_abf.txt") 
test.g = m.sym[i.flow]
tg.a.nb = report.g.cor.invT(test.g,epx.met,pred.a.nb.mag.x,pred.a.nb.spec.x,pred.a.nb.score.g,met.mns,"z_allmeta_corinvT_g_a_nb.txt") 
test.g = m.sym[i.foc]
tg.b.nb = report.g.cor.invT(test.g,epx.met,pred.b.nb.mag.x,pred.b.nb.spec.x,pred.b.nb.score.g,met.mns,"z_allmeta_corinvT_g_b_nb.txt") 
tg.b.f = report.g.cor.invT(test.g,epx.met,pred.b.f.mag.x,pred.b.f.spec.x,pred.b.f.score.g,met.mns,"z_allmeta_corinvT_g_b_f.txt") 





# ===============================================================================      # ** ------ 4 per gene plotting corr - 1/Surv(T) ---- **
test.g = 'CBX5'
gind = match(test.g,rownames(epx.met))
ind.ex = rep(T,length(met.mns))					# * ---- remove surv = 0  ( 5 subjects #828 --> 823)
ind.ex[met.mns==0] = F
gc.mns = met.mns[ind.ex]
gc.st = met.st[ind.ex]
gc.epx = epx.met[gind,][ind.ex]
gc.mag = pred.ab.nb.mag.x[gind,][ind.ex]
gc.spec = pred.ab.nb.spec.x[gind,][ind.ex]
gc.score = pred.ab.nb.score.g[gind,][ind.ex]

# =========================================================================================== * plot *
par(mfrow=c(4,5))

# ---**** Survival **** ---			(1st)
factor.c = 'Exprs'
plot.surv(epx.met[gind,],met.mns,met.st,sprintf("Meta-free Surv: %s-%s",test.g,factor.c))

# ---**** histogram **** ---		(1st)
MFS = gc.mns
x = gc.epx
ind.q = order(x)
ind.st = gc.st[ind.q]
spcol.2 <- colorRampPalette(c("yellow", "orange"))(length(gc.spec))						
spcol.2[ind.st==1] = 'red'
p.sort = summary(x)														  
Q = rep(0,length(x))
Q[x<=p.sort[[2]]] = 1
Q[x>p.sort[[2]]] = 2
Q[x>p.sort[[3]]] = 3
Q[x>p.sort[[5]]] = 4
hist(MFS[Q==1],col='darkgrey',breaks=30,main = sprintf("%s - Q1",factor.c),xlab = 'months')
hist(MFS[Q==2],col=2,breaks=30,main = sprintf("%s - Q2",factor.c),xlab='months')
hist(MFS[Q==3],col=3,breaks=30,main = sprintf("%s - Q3",factor.c),xlab='months')
hist(MFS[Q==4],col=4,breaks=30,main = sprintf("%s - Q4",factor.c),xlab='months')

# ---**** Survival **** ---			(2nd) ******
factor.c = 'Spec'
plot.surv(pred.ab.nb.spec.x[gind,],met.mns,met.st,sprintf("Meta-free Surv: %s-%s",test.g,factor.c))

# ---**** histogram **** ---		(2nd) ******
MFS = gc.mns
x = gc.spec
ind.q = order(x)
ind.st = gc.st[ind.q]
spcol.2 <- colorRampPalette(c("yellow", "orange"))(length(gc.spec))						
spcol.2[ind.st==1] = 'red'
p.sort = summary(x)														  
Q = rep(0,length(x))
Q[x<=p.sort[[2]]] = 1
Q[x>p.sort[[2]]] = 2
Q[x>p.sort[[3]]] = 3
Q[x>p.sort[[5]]] = 4
hist(MFS[Q==1],col='darkgrey',breaks=30,main = sprintf("%s - Q1",factor.c),xlab = 'months')
hist(MFS[Q==2],col=2,breaks=30,main = sprintf("%s - Q2",factor.c),xlab='months')
hist(MFS[Q==3],col=3,breaks=30,main = sprintf("%s - Q3",factor.c),xlab='months')
hist(MFS[Q==4],col=4,breaks=30,main = sprintf("%s - Q4",factor.c),xlab='months')

# ---**** Survival **** ---			(3rd) ******
factor.c = 'SNR'
plot.surv(pred.ab.nb.score.g[gind,],met.mns,met.st,sprintf("Meta-free Surv: %s-%s",test.g,factor.c))

# ---**** histogram **** ---		(3rd) ******
MFS = gc.mns
x = gc.score
ind.q = order(x)
ind.st = gc.st[ind.q]
spcol.2 <- colorRampPalette(c("yellow", "orange"))(length(gc.spec))						
spcol.2[ind.st==1] = 'red'
p.sort = summary(x)														  
Q = rep(0,length(x))
Q[x<=p.sort[[2]]] = 1
Q[x>p.sort[[2]]] = 2
Q[x>p.sort[[3]]] = 3
Q[x>p.sort[[5]]] = 4
hist(MFS[Q==1],col='darkgrey',breaks=30,main = sprintf("%s - Q1",factor.c),xlab = 'months')
hist(MFS[Q==2],col=2,breaks=30,main = sprintf("%s - Q2",factor.c),xlab='months')
hist(MFS[Q==3],col=3,breaks=30,main = sprintf("%s - Q3",factor.c),xlab='months')
hist(MFS[Q==4],col=4,breaks=30,main = sprintf("%s - Q4",factor.c),xlab='months')


# * ------------ **
# * GLM 
# * ------------ **
factor = 'Exprs'
x = gc.epx
xlabs = sprintf("sorted by %s",factor)

plot(1/(gc.mns)~x,col=spcol.2,ylab="1 / Meta-free Surv (months)",xlab= xlabs,main=sprintf("GLM: %s-%s vs MFS",test.g,factor))
abline(glm(1/(gc.mns)~x,family=gaussian),lty=2,lwd=1.5,col='black')
stat = summary(glm(1/(gc.mns)~x,family=gaussian),lty=2,lwd=1.5,col='black')
txt.esti = sprintf("Estimate (GLM) = %.3f",stat[[12]][2])
txt.pv = sprintf("P-value = %.3f",stat[[12]][8])
cat(txt.pv)
legend(min(x),max(1/gc.mns),c("meta-free","meta",txt.esti,txt.pv),cex=0.7,col=c('orange','red'),pch=c(1,1,-1,-1),text.col ='black', bg ='white',bty="n")
