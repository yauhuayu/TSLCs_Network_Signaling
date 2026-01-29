# ----------------------------------------------------------------------------------------------------------------------------------------- #
# * surv lvcom
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

file = c('msk1','msk2','beer4573','shedden-dfci','shedden-hlm','shedden-ummi','bild3141')
dt.files = vector('list',7)

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
for (n in 3:6) {
	cat("Reading network expression data ......\t",file[n],"\n\n")
	
	inf = sprintf("lvcom_nwEpx_%s.txt",file[n])
	data <- read.table(inf, sep="\t",header=T,strip.white=T)
	
	dt.files[[n]]$sample = as.character(data[,1][drop=T])
	dt.files[[n]]$surv.st = as.numeric(data[,2][drop=T])
	dt.files[[n]]$surv.mns = as.numeric(data[,3][drop=T])
	dt.files[[n]]$data = data.matrix(data[,4:ncol(data)])
}
# --------------------------------------------------------------- # 
n=7
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
# *** Surv *** #
epx.surv = e.sym	
surv.st = c()
surv.mns = c()
for (e in 1:6) {
	dd = dt.files[[e]]$data
	dd.1 = dd[,match(e.sym,colnames(dd),nomatch=0)]
	#dd.2 = cbind(dd.1,rep(0,nrow(dd.1)))
	dd.2 = dd.1
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

# ** wt.exprs **
pred.ab.nb.wtEpx.x = epx.surv*wt.ab.nb
pred.ab.nb.wtEpx = apply(pred.ab.nb.wtEpx.x,2,sum)/nrow(epx.surv)
pred.ab.f.wtEpx.x = epx.surv*wt.ab.f
pred.ab.f.wtEpx = apply(pred.ab.f.wtEpx.x,2,sum)/nrow(epx.surv)

pred.a.nb.wtEpx.x = epx.surv*wt.a.nb
pred.a.nb.wtEpx = apply(pred.a.nb.wtEpx.x,2,sum)/nrow(epx.surv)
pred.b.nb.wtEpx.x = epx.surv*wt.b.nb
pred.b.nb.wtEpx = apply(pred.b.nb.wtEpx.x,2,sum)/nrow(epx.surv)
pred.b.f.wtEpx.x = epx.surv*wt.b.f
pred.b.f.wtEpx = apply(pred.b.f.wtEpx.x,2,sum)/nrow(epx.surv)

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


# ================================================================================  				# ** ------ 1 test survival --- ** #
test.ab.nb = report.nw.surv(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score,surv.mns,surv.st) 
test.ab.f = report.nw.surv(pred.ab.epx,pred.ab.f.mag,pred.ab.f.spec,pred.ab.f.score,surv.mns,surv.st) 
test.a.nb = report.nw.surv(pred.a.epx,pred.a.nb.mag,pred.a.nb.spec,pred.a.nb.score,surv.mns,surv.st) 
test.b.nb = report.nw.surv(pred.b.epx,pred.b.nb.mag,pred.b.nb.spec,pred.b.nb.score,surv.mns,surv.st) 
test.b.f = report.nw.surv(pred.b.epx,pred.b.f.mag,pred.b.f.spec,pred.b.f.score,surv.mns,surv.st) 

test = c('epx','mag','spec','snr')
test = data.matrix(cbind(test,test.ab.nb,test.ab.f,test.a.nb,test.b.nb,test.b.f))
write.table(test,"lvcom_surv_nw.txt",sep="\t",col.names =T, row.names=F, quote=F)

# == use wt.Exprs ==
test.ab.nb.2 = report.nw.surv(pred.ab.epx,pred.ab.nb.wtEpx,pred.ab.nb.spec,pred.ab.nb.score,surv.mns,surv.st)
test.ab.f.2 = report.nw.surv(pred.ab.epx,pred.ab.f.wtEpx,pred.ab.f.spec,pred.ab.f.score,surv.mns,surv.st) 
test.a.nb.2 = report.nw.surv(pred.a.epx,pred.a.nb.wtEpx,pred.a.nb.spec,pred.a.nb.score,surv.mns,surv.st)
test.b.nb.2 = report.nw.surv(pred.b.epx,pred.b.nb.wtEpx,pred.b.nb.spec,pred.b.nb.score,surv.mns,surv.st)
test.b.f.2 = report.nw.surv(pred.b.epx,pred.b.f.wtEpx,pred.b.f.spec,pred.b.f.score,surv.mns,surv.st)

test.2 = c('epx','wtEpx','spec','snr')
test = data.matrix(cbind(test.2,test.ab.nb.2,test.ab.f.2,test.a.nb.2,test.b.nb.2,test.b.f.2))
write.table(test,"lvcom_surv_nw_2.txt",sep="\t",col.names =T, row.names=F, quote=F)


# == *** plotting *** ==																			*** ------ all genes (nw) plotting --- **
plot.surv(pred.a.nb.spec,surv.mns,surv.st,'Merged overall survivals by averaged exprs - intermodular hubs')


plot.dyn.nw(pred.b.epx,pred.b.f.mag,pred.b.f.spec,pred.b.f.score,surv.mns,surv.st) 


# == * corr - 1/Surv(T) * ==	------------------------------------------------------------------	* corr - 1/Surv(T) * =============

# ********** none was found significant in linear association with 1/Surv(T) ***********************


# ================================================================================  	# ** ------ 2 per gene all data test survival --- **
test.g = rownames(epx.surv)
tg.ab.nb = report.g.surv(test.g,epx.surv,pred.ab.nb.mag.x,pred.ab.nb.spec.x,pred.ab.nb.score.g,surv.mns,surv.st,"z_all_surv_g_abnb.txt") 
tg.ab.f = report.g.surv(test.g,epx.surv,pred.ab.f.mag.x,pred.ab.f.spec.x,pred.ab.f.score.g,surv.mns,surv.st,"z_all_surv_g_abf.txt") 
test.g = m.sym[i.flow]
tg.a.nb = report.g.surv(test.g,epx.surv,pred.a.nb.mag.x,pred.a.nb.spec.x,pred.a.nb.score.g,surv.mns,surv.st,"z_all_surv_g_a_nb.txt") 
test.g = m.sym[i.foc]
tg.b.nb = report.g.surv(test.g,epx.surv,pred.b.nb.mag.x,pred.b.nb.spec.x,pred.b.nb.score.g,surv.mns,surv.st,"z_all_surv_g_b_nb.txt") 
tg.b.f = report.g.surv(test.g,epx.surv,pred.b.f.mag.x,pred.b.f.spec.x,pred.b.f.score.g,surv.mns,surv.st,"z_all_surv_g_b_f.txt") 



# ================================================================================  	# ** ----- 3 per gene all data cor.test 1 / T(surv) --- **
test.g = rownames(epx.surv)
tg.ab.nb = report.g.cor.invT(test.g,epx.surv,pred.ab.nb.mag.x,pred.ab.nb.spec.x,pred.ab.nb.score.g,surv.mns,"z_all_corinvT_g_abnb.txt") 
tg.ab.f = report.g.cor.invT(test.g,epx.surv,pred.ab.f.mag.x,pred.ab.f.spec.x,pred.ab.f.score.g,surv.mns,"z_all_corinvT_g_abf.txt") 
test.g = m.sym[i.flow]
tg.a.nb = report.g.cor.invT(test.g,epx.surv,pred.a.nb.mag.x,pred.a.nb.spec.x,pred.a.nb.score.g,surv.mns,"z_all_corinvT_g_a_nb.txt") 
test.g = m.sym[i.foc]
tg.b.nb = report.g.cor.invT(test.g,epx.surv,pred.b.nb.mag.x,pred.b.nb.spec.x,pred.b.nb.score.g,surv.mns,"z_all_corinvT_g_b_nb.txt") 
tg.b.f = report.g.cor.invT(test.g,epx.surv,pred.b.f.mag.x,pred.b.f.spec.x,pred.b.f.score.g,surv.mns,"z_all_corinvT_g_b_f.txt") 





# ===============================================================================      # ** ------ 4 per gene plotting corr - 1/Surv(T) ---- **
test.g = 'SHC1'
factor.c = 'Exprs'
gind = match(test.g,rownames(epx.surv))
ind.ex = rep(T,length(surv.mns))					# * ---- remove surv = 0  ( 5 subjects #828 --> 823)
ind.ex[surv.mns==0] = F
gc.mns = surv.mns[ind.ex]
gc.st = surv.st[ind.ex]
gc.epx = epx.surv[gind,][ind.ex]
gc.mag = pred.ab.nb.mag.x[gind,][ind.ex]
gc.spec = pred.ab.nb.spec.x[gind,][ind.ex]
gc.score = pred.ab.nb.score.g[gind,][ind.ex]

