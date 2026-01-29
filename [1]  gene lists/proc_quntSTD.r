# ----------------------------
# date: 2009-04-01
# colon (11 chips)
# ----------------------------

library('affyQCReport')
library('limma')
library('Mfuzz')
library('fdrtool')
source("sub_func.r")


# ----- annotation ----- # 										# --- annotation --- #
ann.dat <-  read.delim("hgU133p2_short.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(ann.dat), " rows, and ", ncol(ann.dat), "columns\n"))
ann.col<-colnames(ann.dat)

# ----- extract 
ann.probe = as.character(ann.dat[,1][drop = T])
g.name = as.character(ann.dat[,3][drop = T])
g.sym = as.character(ann.dat[,4][drop = T])
g.chr = as.character(ann.dat[,5][drop = T])
g.id = as.character(ann.dat[,6][drop = T])



# ---- read-in CEL files ------------------------------------------------------------- proc for 1st: bk, norm, STD ------------ #
test = justRMA()
pt.name = sampleNames(test)

# [1] "HT29-n1.CEL"       "HT29-n2.CEL"       "HT29-n3.cel"       "HT29-n4.cel"       "HT29-n5.cel"       "HT29-sf1.CEL"      "HT29-sf2.CEL"     
# [8] "SW480-n1.CEL"      "SW480-n2.CEL"      "SW480-snail-1.CEL" "SW480-snail-2.CEL"

# ---- pre-process data ---- #
eset.st = standardise(test)
exprs.st = exprs(eset.st)
qt.c = exprs.st[,1]
length(qt.c[is.na(qt.c)])		# *  probes are NAs
exprs.st[is.na(qt.c),]=0

# ** plot
brewer.cols = brewer.pal(8,'Set1')
boxplot(eset.st,col = brewer.cols, ylab = 'Processed log(base2) scale Probe Intensities', xlab = 'Array Names')
hist(eset.st,col = brewer.cols, lty =1, xlab = 'Log(base2) Intensities', lwd = 2)
legend(0,0.2,pt.name,lty=1,col=brewer.cols, lwd=2, cex = 0.7)

# =========== 
probe = rownames(exprs.st)
gr = c(rep('p',5),rep('csc',2),rep('p',2),rep('csc',2))

# * calculate fold changes 
ind.ctrl = grep('p',gr)
ind.tx = grep('csc',gr)
exprs.ctrl = apply(exprs.st[,ind.ctrl],1,mean)
exprs.tx = apply(exprs.st[,ind.tx],1,mean)
fold.df.p = exprs.tx - exprs.ctrl
summary(fold.df.p)



# ------------------------------------ Analysis: 'gene list / modules'  ----------------------- #
# * reading gene lists 
t.dat <-  read.delim("g_probe.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(t.dat), " rows, and ", ncol(t.dat), "columns\n"))
t.probe = as.character(t.dat[,1][drop = T])
t.gsym = as.character(t.dat[,4][drop = T])

index = match(t.probe,ann.probe)

t.exprs = exprs.st[index,c(ind.ctrl,ind.tx)]				# ***
ann.t = ann.dat[index,]
# rownames(t.exprs) = t.gsym

t.fdf = fold.df.p[index]	# fold diff of SAS
summary(t.fdf)


# ---------------- Analysis: 'limma package' -------------- #
samp.gr = factor(c(rep('p',length(ind.ctrl)),rep('csc',length(ind.tx))))
design = model.matrix(~samp.gr)

ck = t.exprs[,1]
t2.fdf = t.fdf[!is.na(ck)]				# exclude  probes
t2.exprs = t.exprs[!is.na(ck),]
t2.probe = t.probe[!is.na(ck)]
t2.gsym = t.gsym[!is.na(ck)]


# * use -- t.exprs
ck.fd = quantile(abs(t2.fdf))
#fdck = 0.5*(ck.fd[4]+ck.fd[5])
fdck = sort(abs(t2.fdf),decreasing=T)[150]
t2.exprs = t2.exprs[abs(t2.fdf)>fdck,]
t2.probe = t2.probe[abs(t2.fdf)>fdck]
t2.gsym = t2.gsym[abs(t2.fdf)>fdck]
t2.fdf = t2.fdf[abs(t2.fdf)>fdck]
t.data = new('ExpressionSet',exprs = t2.exprs)

fit = lmFit(t.data,design)
eb = eBayes(fit)	#  using moderated t-statistic
qqt(eb$t,df=eb$df.prior+eb$df.residual,main="Moderated t")
abline(0,1)						#  * Points off the line may be differentially expressed

eb.p = eb$p.value[,2]
length(eb.p[eb.p<0.05])		# probes:  for p<0.05 ;; 
eb.probe = t2.probe[eb.p<0.05]
ann.eb = ann.dat[match(eb.probe,ann.probe),]



# ------------------ fdrtool -------------------------------------------------------------------------- #
f.probe = t2.probe[eb.p<0.09]
eb.p2 = eb.p[eb.p<0.09]
x = fdrtool(eb.p2,statistic='pvalue',cutoff.method='pct0',pct0=0.3)	# 87 (P<0.05) / 328 (p<0.08) = 0.26

sig.probe = f.probe[x$lfdr<0.05]
q = eb.probe[match(sig.probe,eb.probe,nomatch=0)]			# all of the 87 probes were included
ann.eb2 = ann.dat[match(sig.probe,ann.probe),]
dim(ann.eb2)
# ----------------------------------------------------------------------------------------------------- #



# * heatmap * #
i.mat = t2.exprs[match(sig.probe,eb.probe),]
i.pv = eb.p[match(sig.probe,eb.probe)]
#rownames(i.mat) = eb.probe
rownames(i.mat) = ann.eb2[,4]
i.fdf = t2.fdf[eb.p<0.05][match(sig.probe,eb.probe)]

hmcol = colorRampPalette(c("blue", "white", "red"),space="Lab")(256)
spcol = c(rep('yellow',length(ind.ctrl)),rep('green',length(ind.tx)))
hv = heatmap(i.mat,col=hmcol,ColSideColors = spcol, margins= c(10,15),cexRow = 0.173 + 0.9/log10(nrow(i.mat)))



# * output ordered matrix and genes * #
rlab = hv$rowInd
clab = hv$colInd
order.probe = eb.probe[rlab]
order.mat = i.mat[rlab,clab]
ann.order = ann.eb[rlab,]
order.fdf = i.fdf[rlab]
order.pv = i.pv[rlab]

out = cbind(order.probe,order.fdf,order.pv,ann.order,order.mat)
write.table(out,'out_cn_20090401.txt',sep = '\t',quote=F,col.names =T ,row.names=F)


