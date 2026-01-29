# ----------------------------------------------------------------------------------------------------------------------------------------- #
# msk1 lvcom network
# * surv
# ------------------------------------------------------------------------------------------------------------------------------------------
library('survival')
library('affyQCReport')
library('limma')
library('Mfuzz')
library(labdsv)
source("sub_func.r")
source("method.r")
source("method_surv.r")
source('proc_g_plot.r')

# ----- probe-CA-specific ----- # 																					# --- annotation --- #
mark <-  read.delim("lvcom_cnvt.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(mark), " rows, and ", ncol(mark), "columns\n"))

m.sym = as.character(mark[,1][drop=T])
m.lvpb = as.character(mark[,7][drop=T])
m.n.type = as.numeric(mark[,2][drop=T])
m.n.nb = as.numeric(mark[,3][drop=T])
m.n.foc = as.numeric(mark[,4][drop=T])
m.a.pbs = as.character(mark[,9][drop=T])


# -------------------- read-in surv exprs ------------------------------- # 			
# msk1.txt
# ----------------------------------------------------------------------- #			
data = read.table("msk1_gerald.txt",sep="\t",header=T,strip.white=T)
cat("reading file: \n")
cat("check data:	",dim(data),"\n")
pt = colnames(data)
pt2 = gsub("\\.CEL","",pt)
pt2 = gsub("_Rep","",pt2)
pt2 = gsub("NCI_U133A_","",pt2)

class = read.delim("msk1_sdrf.txt",sep="\t",header=T,strip.white=T)
dim(class)
sample = as.character(class[,1][drop=T])
met.st = as.numeric(class[,2][drop=T])
surv.st = as.character(class[,3][drop=T])
surv.st[surv.st=="A"] = 0
surv.st[surv.st=="D"] = 1
surv.st = as.numeric(surv.st)
met.mns = as.numeric(class[,4][drop=T])
surv.mns = as.numeric(class[,5][drop=T])



# ---- transform data type 'list' to 'data.matrix' ----- 
s.exprs = data
colnames(s.exprs) = pt2
s.pb = rownames(s.exprs)
k = dim(s.exprs)
exprs = data.matrix((unlist(s.exprs)))
exprs = as.numeric(exprs)
dim(exprs) = k
dim(exprs)

exprs.n = normalizeQuantiles(exprs)
n.set = new('ExpressionSet',exprs = exprs.n)
n.set = standardise(n.set)
exprs.st = exprs(n.set)
qt.c = exprs.st[,1]
length(qt.c[is.na(qt.c)])		# *  probes are NAs
exprs.st[is.na(qt.c),]=0
colnames(exprs.st) = pt2
rownames(exprs.st) = s.pb

# ----------------------------------- #
# match p2 to u133A probes
# ----------------------------------- #
cnvt.ann = mark[!is.na(m.a.pbs),]
dim(cnvt.ann)
cnvt.pbs = m.a.pbs[!is.na(m.a.pbs)]
cnvt.n.type = m.n.type[!is.na(m.a.pbs)]
cnvt.n.nb = m.n.nb[!is.na(m.a.pbs)]
cnvt.n.foc = m.n.foc[!is.na(m.a.pbs)]
cnvt.gsym = m.sym[!is.na(m.a.pbs)]


# ******************** extract converted pbs Exprs ******************* #
cnvt.exprs = array(9999,dim = c(nrow(cnvt.ann),ncol(exprs.st)))
for (w in 1:nrow(cnvt.ann)) {
	w.pbs = unlist(strsplit(cnvt.pbs[w],","))
	w.pbs = trimWhiteSpace(w.pbs)
	cat("ck\t",w,"\t",w.pbs,"\n")
	if(length(w.pbs)==1) {
		cnvt.exprs[w,] = exprs.st[match(w.pbs,s.pb),]
	}
	if (length(w.pbs)>1) {
		w.epx = exprs.st[match(w.pbs,s.pb),]
		w.mx = apply(abs(w.epx),2,max)
		cnvt.exprs[w,] = w.mx
	}
}
colnames(cnvt.exprs) = pt2
pt2[match(sample,pt2)][1:10]


# ** ---- re-order the exprs file according to the surv sdrf file ---- ** #
cnvt.epx = cnvt.exprs[,match(sample,pt2)]			# * ind.dd and alv according to sdrf
rownames(cnvt.epx) = cnvt.gsym
# ************************* #
ind.alv = grep(0,surv.st)	
ind.dd = grep(1,surv.st)
#ind.alv = grep(0,met.st)	
#ind.dd = grep(1,met.st)
s.exprs.v = cnvt.epx[,c(ind.alv,ind.dd)]
dim(s.exprs.v)

# ************************* #
surv.mns.x = surv.mns[c(ind.alv,ind.dd)]
surv.st.x = surv.st[c(ind.alv,ind.dd)]
#surv.mns.x = met.mns[c(ind.alv,ind.dd)]
#surv.st.x = met.st[c(ind.alv,ind.dd)]
cbind(c(ind.alv,ind.dd),surv.st.x,surv.mns.x,sample[c(ind.alv,ind.dd)],colnames(s.exprs.v))[1:4,]


# ****** output files ****** #
length(surv.st.x[surv.st.x==1])
t.s.epx = t(s.exprs.v)
out.mat = cbind(colnames(s.exprs.v),surv.st[c(ind.alv,ind.dd)],surv.mns[c(ind.alv,ind.dd)],met.st[c(ind.alv,ind.dd)],met.mns[c(ind.alv,ind.dd)],t.s.epx)
write.table(out.mat,"lvcom_nwEpx_msk1.txt",sep="\t", col.names =T, row.names=F, quote=F)



# ================================================================================================================================== #
# test survival using consensus network: DNMT3A / SIN3A
# ================================================================================================================================== #

i.foc = grep(1,m.n.type)
i.flow = grep(2,m.n.type)
e.sym = m.sym[c(i.foc,i.flow)]

s.exprs.v = cnvt.epx[,c(ind.alv,ind.dd)]
epx.surv = s.exprs.v[match(e.sym,rownames(s.exprs.v),nomatch=0),]


# *** modify variable ***
surv.mns = met.mns[c(ind.alv,ind.dd)]
surv.st = met.st[c(ind.alv,ind.dd)]

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
# **** Global statistics **** # 
#w = cbind(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score)
#w = cbind(pred.ab.epx,pred.ab.f.mag,pred.ab.f.spec,pred.ab.f.score)
#w = cbind(pred.a.epx,pred.a.nb.mag,pred.a.nb.spec,pred.a.nb.score)
#w = cbind(pred.b.epx,pred.b.nb.mag,pred.b.nb.spec,pred.b.nb.score)
#w = cbind(pred.b.epx,pred.b.f.mag,pred.b.f.spec,pred.b.f.score)
#r.test.out.1 = report.nw.surv(pred.ab.epx,pred.ab.nb.mag,pred.ab.nb.spec,pred.ab.nb.score,surv.mns,surv.st) 

# ------------------------------------------------------------------------------------------------------------------------------------------
library(geneplotter)
for (gtest in 1:length(e.sym)) {
	#test.g = 'SIN3A'
	test.g = m.sym[gtest]
	gind = match(test.g,rownames(epx.surv))

# **** =========== ab.nb =========== ****
gc.epx = epx.surv[gind,]
gc.mag = pred.ab.f.mag.x[gind,]
gc.spec = pred.ab.f.spec.x[gind,]
gc.score = pred.ab.f.score.g[gind,]

# **** statistics **** # 
#w = cbind(gc.epx,gc.mag,gc.spec,gc.score)
#w0 = w[surv.st==0,]
#w0.mns = surv.mns[surv.st==0]
#w1 = w[surv.st==1,]
#w1.mns = surv.mns[surv.st==1]
#type = c('Exprs','Magnitude','Spectrum','SNR')
cat("Gene:  ",test.g,"\t\t")
xx = report.nw.surv(gc.epx,gc.mag,gc.spec,gc.score,surv.mns,surv.st) 
cat(xx,"\n")
} #=== end of for gtest loop

# * --------------------
# * Survival analysis
# * Histogram
# * --------------------
par(mfrow=c(2,2))
for (y in 1:4) {
xtitle = sprintf("G-%s-%s",test.g,type[y])
plot.surv(w[,1],surv.mns,surv.st,xtitle)
}
