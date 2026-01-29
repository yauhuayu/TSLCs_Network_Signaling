# ----------------------------------------------------------------------------------------------------------------------------------------- #
# * surv msk 1 - gerald
# * use lng_cnvt.txt 
# ------------------------------------------------------------------------------------------------------------------------------------------
library('survival')
library('affyQCReport')
library('limma')
library('Mfuzz')
library('fdrtool')
library(labdsv)
source("sub_func.r")
source("method.r")

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

exprs.n = normalize.quantiles(exprs)
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
cbind(c(ind.alv,ind.dd),surv.st.x,surv.mns.x,sample[c(ind.alv,ind.dd)])[1:4,]


# ****** output files ****** #
length(surv.st.x[surv.st.x==1])
t.s.epx = t(s.exprs.v)
out.mat = cbind(colnames(s.exprs.v),surv.st[c(ind.alv,ind.dd)],surv.mns[c(ind.alv,ind.dd)],met.st[c(ind.alv,ind.dd)],met.mns[c(ind.alv,ind.dd)],t.s.epx)
write.table(out.mat,"lng_nwEpx_msk1.txt",sep="\t", col.names =T, row.names=F, quote=F)





# ================================================================================================================================== #
# calculate prediction of survival 
# by topology																	  
#
# ================================================================================================================================== #
i.foc = grep(1,cnvt.n.type)
i.flow = grep(2,cnvt.n.type)

s.exprs.v = cnvt.epx[,c(ind.alv,ind.dd)]
s.exprs.v = s.exprs.v[c(i.foc,i.flow),]

#  ========== ALL hub GENES ============ #                                       
wt.fl = rep(0,nrow(s.exprs.v))
#wt.fl[i.flow] = cnvt.n.nb[i.flow]
#wt.fl[i.foc] = cnvt.n.nb[i.foc]
wt.fl[1:length(i.foc)] = cnvt.n.foc[i.foc]
wt.fl[(length(i.foc)+1):length(wt.fl)] = cnvt.n.nb[i.flow]

# --- eg for gene specific 
eg.spec.1= apply(abs(s.exprs.v),2,k.b.wt.dis.G,wt.fl)/(nrow(s.exprs.v)-1)
eg.spec.2= apply(abs(s.exprs.v),2,k.b.wt.dis.G.contra,wt.fl)/(nrow(s.exprs.v)-1)
eg.spec = 0.5*(eg.spec.1 + eg.spec.2)

fl.spec = apply(eg.spec,2,sum)/nrow(s.exprs.v)
fl.spec.ck = apply(eg.spec.1,2,sum)/nrow(s.exprs.v)
#cbind(fl.spec.ck,fl.spec)									# the original eg.spec or fl.spec is the same when try to sum symmetric kb.dist

eg.mag = abs(s.exprs.v)*wt.fl
fl.mag = apply(eg.mag,2,sum)/nrow(s.exprs.v)
fl.score = fl.mag/fl.spec
eg.score = eg.mag/eg.spec

# ******
wt.epx = s.exprs.v*wt.fl
t.epx = apply(wt.epx,2,sum)/nrow(s.exprs.v)					# without abs with topological info		

wt.t1 = wt.fl
wt.t1[wt.fl!=0] = 1
t1.epx = apply(s.exprs.v*wt.t1,2,sum)/nrow(s.exprs.v)		# sum / ave of expression values only (no topological information)		

length(wt.t1[wt.t1==1])
length(surv.st.x [surv.st.x ==1])
length(met.st[met.st==1])
length(surv.st[surv.st==1])



# ------------------------- test survival ----------------------------------------------------------------------------- Survival ---------- # 			
# generate surv plot																
# ----------------------------------------------------------------------------------------------------------------------------------------- #	
out.score = t1.epx
p.sort = summary(out.score)														    				# ** --- 1 --- ** #

p.gen = rep(0,length(out.score))
p.gen[out.score<=p.sort[[2]]] = 1
p.gen[out.score>p.sort[[2]]] = 2
p.gen[out.score>p.sort[[3]]] = 3
p.gen[out.score>p.sort[[5]]] = 4

svt.data = as.data.frame(cbind(surv.mns.x,surv.st.x,p.gen))
dim(svt.data)
svt.data = as.data.frame(as.list(svt.data))
dim(svt.data)

fit <- survfit(Surv(surv.mns.x,surv.st.x)~  p.gen, data = svt.data, type = 'kaplan-meier')
plot(fit, col = 1:4,xlab='years',ylab='survival')
title(main='Survival curves grouped by scores derived from hubs',cex.main=1)
legend(60, 0.95, c("q1","q2","q3","q4","p-value = 0.22"), cex=0.7, col = 1:4,lty = c(rep(1,4),0), text.col = 1:4, merge = TRUE, bg = 'white')
survdiff(Surv(surv.mns.x,surv.st.x)~  p.gen, data = svt.data)



# ================================================================================ 					 # ** --- 2 --- ** #
p.gen2 = rep(0,length(out.score))
p.gen2[p.gen==3] = 1
p.gen2[p.gen==4] = 1

svt.data = as.data.frame(cbind(surv.mns.x,surv.st.x,p.gen2))
dim(svt.data)
svt.data = as.data.frame(as.list(svt.data))
dim(svt.data)

fit <- survfit(Surv(surv.mns.x,surv.st.x)~  p.gen2, data = svt.data, type = 'kaplan-meier')
plot(fit,lty = 2:3,mark.time=F)
points(fit[1],pch = 3,col=2)
points(fit[2],pch = 4,col=4)
title(main='Survival curves grouped by scores derived from the intra-modular hubs',cex.main=1)
legend(70, 0.95,c("q(1,3,4)","q(2)","p-value = 0.03"),cex=0.7,col=c(2,4,1),lty =c(2,3,0),pch=c(3,4,-1),text.col = 'black', merge = TRUE, bg ='white')
survdiff(Surv(surv.mns.x,surv.st.x)~  p.gen2, data = svt.data)

