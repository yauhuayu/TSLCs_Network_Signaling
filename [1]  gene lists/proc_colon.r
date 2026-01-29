# ----------------------------------------------------------------------------------------------------------------------------------------------- #
# *  note: colon (11 chips)
# ----------------------------------------------------------------------------------------------------------------------------------------------- #

library('affyQCReport')
library('limma')
library('Mfuzz')
library('fdrtool')
source("sub_func.r")
source("method.r")

# ----- annotation ----- # 										# --- annotation --- #
ann.dat <-  read.delim("hgu133p2_ann.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(ann.dat), " rows, and ", ncol(ann.dat), "columns\n"))
ann.col<-colnames(ann.dat)
# "ann.probe" "g.id.num"  "out.gid"   "g.sym"     "g.chr"     "g.omim"    "n.ann.tsp" "g.name"

# ----- extract 
a.pb = as.character(ann.dat[,1][drop = T])
a.gid.n = as.numeric(ann.dat[,2][drop=T])
a.gids = as.character(ann.dat[,3][drop = T])
a.sym = as.character(ann.dat[,4][drop = T])
a.chr = as.character(ann.dat[,5][drop = T])
a.omim = as.character(ann.dat[,6][drop = T])
a.gtsp.n = as.numeric(ann.dat[,7][drop=T])
a.gname = as.character(ann.dat[,8][drop = T])



# ---- read-in CEL files ------------------------------------------------------------- proc for 1st: bkgn, norm, STD ------------ #
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
exprs.tx[exprs.tx==0] = 0.0001									# So that, the var.tx is not NA, and of large values
sd.tx = apply(exprs.st[,ind.tx],1,sd)					
var.tx = sd.tx/abs(exprs.tx)
length(var.tx[is.na(var.tx)])		
length(var.tx[var.tx==0])										# only 2 CSC samples, 2682 probes with same values, var.tx == 0

fold.df.p = exprs.tx - exprs.ctrl
summary(fold.df.p)



# ------------------------------------ Process: 'gene list / modules'  ----------------------- #
# * reading gene lists 
t.dat <-  read.delim("ca_probe.txt", sep="\t",header=F)
cat(c("The annotation data contains ", nrow(t.dat), " rows, and ", ncol(t.dat), "columns\n"))
ca.pbs = trimWhiteSpace(as.character(t.dat[,1][drop = T]))


# * reading hgu133p2_gid
g2pb = read.delim("hgu133p2_gid.txt", sep="\t",header=F,strip.white=T)
g2pb.gid = trimWhiteSpace(as.character(g2pb[,1][drop=T]))
g2pb.pbs = trimWhiteSpace(as.character(g2pb[,6][drop=T]))


# * reading "lung cancer" - specific alternative splicing genes
br.altsp = read.delim("cn_altsp.txt", sep="\t",header=T,strip.white=T)
br.asp.gid = trimWhiteSpace(as.character(br.altsp[,3][drop=T]))
br.asp.g2 = gsub("\\*","",br.asp.gid)

br.ind = match(br.asp.g2,g2pb.gid,nomatch=0)
length(br.ind[br.ind==0])							# 8 geneIDs not matched to the final annoated "hgu133p2_gid.txt"
br.asp.g2[br.ind==0]								# * some are hypothetical proteins or pseudogenes

br.g = g2pb[br.ind,]
br.tot = sum(br.g[,2])								# 1464
br.asp.pbs = trimWhiteSpace(as.character(br.g[,6][drop=T]))
br.asp.pbs = trimWhiteSpace(unlist(strsplit(br.asp.pbs,",")))
length(br.asp.pbs)

t.probe = removeRepeat(c(ca.pbs,br.asp.pbs))
length(t.probe)										# 13558 (13067 + 491)


# ---------------------------------- Analysis: 'limma - two groups'  -------------------------------- #
index = match(t.probe,probe)

t.exprs = exprs.st[index,c(ind.ctrl,ind.tx)]				# ***
ann.t = ann.dat[index,]
t.gsym = a.sym[index]
# rownames(t.exprs) = t.gsym

t.fdf = fold.df.p[index]	# fold diff of SAS
t.var.tx = var.tx[index]
summary(t.fdf)


# ** design of limma - grouping ** #
samp.gr = factor(c(rep('p',length(ind.ctrl)),rep('csc',length(ind.tx))))
design = model.matrix(~samp.gr)

ck = t.exprs[,1]
t2.fdf = t.fdf[!is.na(ck)]				# exclude  probes
t2.exprs = t.exprs[!is.na(ck),]
t2.probe = t.probe[!is.na(ck)]
t2.gsym = t.gsym[!is.na(ck)]


# * use -- t.exprs
ck.fd = quantile(abs(t2.fdf))
fdck = 0.5*(ck.fd[4]+ck.fd[5])
#fdck = sort(abs(t2.fdf),decreasing=T)[150]
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
ann.eb = ann.dat[match(eb.probe,a.pb),]



# ------------------ fdrtool -------------------------------------------------------------------------- #
f.probe = t2.probe[eb.p<0.09]
eb.p2 = eb.p[eb.p<0.09]
x = fdrtool(eb.p2,statistic='pvalue',cutoff.method='pct0',pct0=0.3)		

sig.probe = f.probe[x$lfdr<0.05]							# diff.exprs probes											# -- A -- #
q = eb.probe[match(sig.probe,eb.probe,nomatch=0)]			
ann.eb2 = ann.dat[match(sig.probe,a.pb),]
dim(ann.eb2)
# ----------------------------------------------------------------------------------------------------- #


# ck
x = read.delim("out_cn_20090501.txt", sep="\t",header=T)
dim(x)
x.pb = as.character(x[,1][drop = T])
x.ind = match(x.pb,sig.probe,nomatch = 0)
length(x.ind[x.ind==0])                       				# the output file on 20090501 was found sig using 0.875% and fdr = 0.05



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
write.table(out,'out_cn_20090522.txt',sep = '\t',quote=F,col.names =T ,row.names=F)




# ********************************************************************************************************************************************* #
#
#												---- build networks -----
# ********************************************************************************************************************************************* #
out.gid = as.character(out[,6][drop=T])
out.tot.gid = sum(as.numeric(out[,5][drop=T]))
out.gid.x = unlist(strsplit(out.gid,"   "))

thrombin = read.delim("thrombin.txt", sep="\t",header=T)
thr.gid = as.character(thrombin[,4][drop=T])

r.gid = removeRepeat(c(thr.gid,out.gid.x))
length(r.gid)

# look into lung-specific !
r.altsp.ind = match(r.gid,lng.asp.g2,nomatch=0)

r.g2pb = g2pb[match(r.gid,g2pb.gid),]
r.sym = as.character(r.g2pb[,3][drop=T])

# ---- interaction network ---------- from gene.Rxn ------------------------ #
g.rxn <-  read.delim("hs_rxn_g.txt", sep="\t",header=T)
cat(c("The interaction data contains ", nrow(g.rxn), " rows, and ", ncol(g.rxn), "columns\n"))
col.grxn = c("id.a","id.b","pmid","txt","db")
colnames(g.rxn) = col.grxn

casi.r2 = ext.list.M2(r.gid,g.rxn)		# *  rxns (nrow)  --> both of the interaction node are in the gene list


# -------------- node.summary -------------- #							
gene.a =as.character(casi.r2[,1][drop=T])
sym.a = get.sym(gene.a,r.gid,r.sym)
gene.b =as.character(casi.r2[,2][drop=T])
sym.b = get.sym(gene.b,r.gid,r.sym)
casi.rxn = cbind(sym.a,sym.b,casi.r2)
ind = rep(T,length(sym.a))
ind[sym.a==sym.b]=F
casi.rxn2 = casi.rxn[ind,]			# * remove self-directing links  (  nodes,  links)

# ***** 
r.g.txt = sprintf("%s :: %s",casi.rxn2[,5],casi.rxn2[,6])
casi.rxn2 = casi.rxn2[,c(1:4,7)]
casi.rxn2 = cbind(casi.rxn2,r.g.txt)

write.table(casi.rxn2,"lung_nw_grxn.txt",sep="\t", col.names =T, row.names=F, quote=F)





# ========================================================================================================================================
node.a = as.character(casi.rxn2[,1][drop=T])
node.b = as.character(casi.rxn2[,2][drop=T])
label = removeRepeat(c(node.a,node.b))
net.m = rxn2mat(node.a,node.b,label)[[1]]

out.nbhood = c()
for (i in 1:length(label)) {
	nb = get.nb(i,net.m,label)
	out.nbhood = append(out.nbhood,length(nb))
	cat(c("ck num: length of out.nbhood= ",length(out.nbhood),"\n"))
}

# ck list
ck.s = out.nbhood		# search Q starting pts
ck.s[out.nbhood==1]=0		# remove peri genes from search Q


out.foc.cp = c()
for (i in 1:length(label)) {
	if (ck.s[i]>0){
		#foc.cp = cut.dist.cp(i,net.m,label,out.nbhood,ck.s,label[i])
		foc.cp = focal.pth.cp(i,net.m,label,out.nbhood,ck.s,label[i])
		out.foc.cp = append(out.foc.cp,foc.cp)
	}
	if (ck.s[i]==0) {
		foc.cp=0
		out.foc.cp = append(out.foc.cp,foc.cp)
	}
	cat(c("ck length of out.foc= ",length(out.foc.cp),"\n"))
}

# ===== ------------ ===== #
# grouping of node										# **
# ===== ------------ ===== #

cat("total number of genes:\t",length(label),"\n")
ind.foc = rep(F,length(label))
ind.foc[out.foc.cp>0]=T
foc.g = label[ind.foc]

ind.flow = rep(T,length(label))
ind.flow[out.foc.cp>0]=F
ind.flow[out.nbhood==1]=F
flow.g = label[ind.flow]

ind.peri = rep(F,length(label))
ind.peri[out.nbhood==1]=T
peri.g = label[ind.peri]

node.type = rep(0,length(label))
node.t.name = rep(NA,length(label))

focality.ind = match(foc.g,label)   		# ck: 	
flow.ind = match(flow.g,label)			# ck: 				
peri.ind = match(peri.g,label)			# ck: 					

node.type[focality.ind] = 1
node.t.name[focality.ind] = "focality"
node.type[flow.ind] = 2
node.t.name[flow.ind] = "flow"
node.type[peri.ind] = 3
node.t.name[peri.ind] = "peri"

node.casi2 = r.g2pb[match(label,r.sym),]
node.sum = cbind(label,node.casi2,node.type,node.t.name,out.nbhood,out.foc.cp)		# only match one line !!

write.table(node.sum,"net_node_hprd.txt",sep="\t", col.names =T, row.names=F, quote=F)


