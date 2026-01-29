# ----------------------------------------------------------------------------------------------------------------------------------------- #
# proc3 generate gene list lv_com
# ----------------------------------------------------------------------------------------------------------------------------------------- #
source("sub_func.r")
library(labdsv)

# ----- annotation ----- # 																							# --- annotation --- #
ann.dat <-  read.delim("gC_cord_FINAL.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(ann.dat), " rows, and ", ncol(ann.dat), "columns\n"))

ann.probe = as.character(ann.dat[,1][drop = T])
g.sym = as.character(ann.dat[,10][drop = T])
g.id = as.character(ann.dat[,9][drop = T])


# ------------------------- read-in pre-processed files --------------------------------------------------------------- read-in ---------- # 			
# quntSTD_ALL.txt
# ---------------------------------------------------------------------------------------------------------------------------------------- #			
data = read.table("quntSTD_ALL_z.txt",sep="\t",header=T,strip.white=T)
cat("reading file: \n")
cat("check data:	",dim(data),"\n")
probe = as.character(data[,1][drop=T])
data = data[,2:ncol(data)]
pt = colnames(data)
pt.name =  gsub(".CEL", "", pt)
pt.name = gsub(".cel", "", pt.name)
pt.name = gsub("\\.", "_", pt.name)

gr = c(rep('csc',4),rep('p',3),'csc','p','csc',rep('p',3),'csc',rep('p',2),
       'csc',rep('Esc',10),rep('p',5),rep('csc',2),'p','csc',rep('p',3),
	   rep('csc',2),rep('Esc',3),rep('csc',4),rep('p',3),rep('Esc',3),
	   rep('p',3),rep('csc',6),rep('Esc',3),rep('p',2),rep('csc',2),
	   rep('Esc',3),'p','csc',rep('Esc',3))


ind.p = grep('p',gr)		#27
length(ind.p)
ind.csc = grep('csc',gr)	#26
length(ind.csc)
ind.Esc = grep('Esc',gr)	#25
length(ind.Esc)

 
# ---- transform data type 'list' to 'data.matrix' ----- 
k = dim(data)
exprs = unlist(data)
dim(exprs) = k
dim(exprs)
colnames(exprs) = pt.name
rownames(exprs) = probe

# --- all genes in the out_cord_FINAL --- #
index = match(ann.probe,probe)			 	#
exprs.v = data.matrix(exprs[index,])
dim(exprs.v)

# ------------------------- test Esc vs p vs csc ---------------------------------------------------------------------- heatmap ---------- # 			
# [1] "v.pb.csc"      "set.a.pb"      "lv_set.a"      "nov_set.a"    
# [5] "lv_csc_nSetA"  "g.lv7.ALL"     "g.lv6.csc.Esc" "g.lv5.p.csc"  
# [9] "g.lv3.csc"     "tcnd"          "cor.mean.csc"  "cor.mean.p"   
# [13] "cor.mean.Esc" 
# ---------------------------------------------------------------------------------------------------------------------------------------- #
ind.m = ann.dat[,c(2,15:17)]
nu = ind.m[,1]
nu[nu==1]=0
nu[nu>1] = 1									# concordance (tis-cond >=2)
i.lv = as.character(ann.dat[,14][drop=T])
i.lv[i.lv=="x"]=1								# extract lv
i.lv[i.lv=="com"]=1								# extract 'com' in gA_tCnd
i.lv[i.lv!="1"] = 0
i.lv = as.numeric(i.lv)
i.lv = i.lv+nu

gsb.p = exprs.v[i.lv==2,ind.p]
gsb.csc = exprs.v[i.lv==2,ind.csc]
gsb.Esc = exprs.v[i.lv==2,ind.Esc]

im.ALL = cbind(gsb.p,gsb.csc,gsb.Esc)
sp.ALL = c(rep('green',length(ind.p)),rep('red',length(ind.csc)),rep('yellow',length(ind.Esc)))
i.ooc = grep("oocyte",colnames(im.ALL))
sp.ALL[i.ooc] = 'black'
cbind(colnames(im.ALL),sp.ALL)

gr.cscEsc = c(rep(1,length(ind.csc)),rep(0,length(ind.Esc)))
gr.pCsc = c(rep(-1,length(ind.p)),rep(1,length(ind.csc)))
gr.ALL = c(rep(-1,length(ind.p)),rep(1,length(ind.csc)),rep(0,length(ind.Esc)))

i.mat = im.ALL
spcol = sp.ALL
i.gr = gr.ALL
ann.lv = ann.dat[i.lv==2,]
rownames(i.mat) = ann.lv[,10]
dim(i.mat)

p.mat = t(i.mat)
p.dist = dist(p.mat)
p.dist[p.dist==0] = 0.00001
#ss = bestnmds(p.dist,k=3,itr = 20,maxit = 100)
ss = nmds(p.dist,k=3)
plot(ss,col = sp.ALL)
title(main='Kruskal\'s Non-metric Multidimensional Scaling')
legend(-7, -7,c("TSLCs","PTCs","hESCs","oocytes"),cex=0.7,col=c('red','green','yellow','black'),pch=1,text.col = 'grey6', bg ='white')


hmcol = colorRampPalette(c("green", "white", "red"),space="Lab")(2048)
hv = heatmap(i.mat,col=hmcol,ColSideColors = spcol, margins= c(10,15),cexRow = 0.173 + 0.9/log10(nrow(i.mat)))

library(marray)
x = as.vector(i.mat)
xx = seq(min(x),max(x),by = 0.01)
maColorBar(hmcol,scale=xx,k=15)

# * output * #
rlab = hv$rowInd
clab = hv$colInd
lv.pb = ann.probe[i.lv==2]
order.probe = lv.pb[rlab]
order.mat = i.mat[rlab,clab]
ann.lv = ann.dat[i.lv==2,]
ann.order = ann.lv[rlab,]
out = cbind(order.probe,ann.order,order.mat)

write.table(out,"lvcom_cord_heatmap.txt",sep = '\t',quote=F,col.names =T ,row.names=F)

