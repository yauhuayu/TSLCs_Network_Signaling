# ----------------------------------------------------------------------------------------------------------------------------------------------- #
# proc4 generate lv_com nw
# ----------------------------------------------------------------------------------------------------------------------------------------------- #

library('limma')
source("sub_func.r")
source("method.r")

# ----- annotation ----- # 																							# --- annotation --- #
a1.dat <-  read.delim("lvcom_cord.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(a1.dat), " rows, and ", ncol(a1.dat), "columns\n"))
a1.pbs = as.character(a1.dat[,1][drop = T])
a1.gname = as.character(a1.dat[,14][drop = T])
a1.gsym= as.character(a1.dat[,11][drop = T])
a1.gid = as.character(a1.dat[,10][drop = T])
a1.gid = unlist(strsplit(a1.gid,"   "))
a1.gid = trimWhiteSpace(a1.gid)

a2.dat <-  read.delim("lvcom_icord.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(a2.dat), " rows, and ", ncol(a2.dat), "columns\n"))
a2.pbs = as.character(a2.dat[,2][drop = T])
a2.gname = as.character(a2.dat[,11][drop = T])
a2.gsym= as.character(a2.dat[,9][drop = T])
a2.gid = as.character(a2.dat[,7][drop = T])

nr.gid = removeRepeat(c(a1.gid,a2.gid))
length(nr.gid)
ov.gid = outRepeat(c(a1.gid,a2.gid))
length(ov.gid)
ind.g = rep("x",length(nr.gid))
ind.g[match(a1.gid,nr.gid,nomatch=0)] = "lvcom"
ind.g[match(a2.gid,nr.gid,nomatch=0)] = "ipa"
ind.g[match(ov.gid,nr.gid)] = "focus"
cbind(nr.gid,ind.g)

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

# ******************************************************************************************************************************************** #
# 
# --- find the max fdf probe to represent "nr.gid" in the network --- #
#												
# ********************************************************************************************************************************************* #

#  * reading hgu133p2_gid *  #
g2pb = read.delim("hgu133p2_gid.txt", sep="\t",header=F,strip.white=T)
g2pb.gid = trimWhiteSpace(as.character(g2pb[,1][drop=T]))
g2pb.sym = trimWhiteSpace(as.character(g2pb[,3][drop=T]))
g2pb.pbs = trimWhiteSpace(as.character(g2pb[,6][drop=T]))

gid.pb = g2pb.pbs[match(nr.gid,g2pb.gid)]
gid.g2pb = g2pb[match(nr.gid,g2pb.gid),]
gid.sym = g2pb.sym[match(nr.gid,g2pb.gid)]

t.probe = c()
t.fdf = c()

for(i in 1:length(nr.gid)) {
    pat = sprintf("^%s$", nr.gid[i])
	s.p = unlist(strsplit(gid.pb[[i]],","))
	s.p = trimWhiteSpace(s.p)
	s.lg = length(s.p)
	cat(pat,"\n",s.p,"\n")
	if (s.lg>1) {
		e.p = apply(exprs[match(s.p,probe),ind.p],1,mean)
		e.csc = apply(exprs[match(s.p,probe),ind.csc],1,mean)
		e.fdf = e.csc-e.p
		ck.fdf = abs(e.fdf)
		t.p = s.p[ck.fdf==max(ck.fdf)]
		t.probe = append(t.probe,t.p)
		t.fdf = append(t.fdf,e.fdf[ck.fdf==max(ck.fdf)])
	}
	if (s.lg ==1) {
		e.p = mean(exprs[match(s.p,probe),ind.p])
		e.csc = mean(exprs[match(s.p,probe),ind.csc])
		e.fdf = e.csc-e.p
		t.probe = append(t.probe,s.p)
		t.fdf = append(t.fdf,e.fdf)
	}
}

t.out1 = cbind(t.probe,t.fdf,nr.gid,ind.g)
t.out2 = g2pb[match(nr.gid,g2pb.gid),]
t.out2 = cbind(t.probe,t.fdf,ind.g,t.out2)

# --------------------------------------------------------------------------------------------------------------------------------------------
# ======= extract exprs by tissue-cond type
# --------------------------------------------------------------------------------------------------------------------------------------------
index = match(t.probe,probe)
exprs.v = data.matrix(exprs[index,])

# ***************** lung *						# 5 vs 3 = 8
i.lng.p = c(5:7)
i.lng.csc = c(1:4,8)
pt.name[i.lng.csc]
# ***************** atrt *						# 5 vs 2 = 7
i.atrt.p = c(9,11:13,35)
i.atrt.csc = c(10,36)
pt.name[i.atrt.csc]
# ***************** oral *						# 5 vs 8 = 13
i.or.p = c(15,16,55:57)
i.or.csc = c(14,17,58:63)
pt.name[i.or.csc]
# ***************** colon *						# 7 vs 4 = 11
i.cn.p = c(28:32,67,68)
i.cn.csc = c(33,34,69,70)
pt.name[i.cn.csc]
# ***************** breast *					# 3 vs 2 = 5
i.br.p = c(37:39)
i.br.csc = c(40,41)
pt.name[i.br.csc]
# ***************** gbm *						# 4 vs 5 = 9
i.gbm.p = c(49:51,74)
i.gbm.csc = c(45:48,75)
pt.name[i.gbm.csc]
# ***************** Snail -- FaDu, MCF7, SW480										# 7 vs 6 = 13
i.snail.p = c(5:7,15,16,37:39,67,68)
i.snail.csc = c(8,14,17,40,41,69,70)
pt.name[i.snail.p]
# **************** Serum-Free Sphere Formation ---- ATRT\HTB, SAS, HT29, U87	  	# 11 vs 11 = 22
i.sf.p = c(9,35,28:32,74)
i.sf.csc = c(10,36,33,34,75)
pt.name[i.sf.p]

#***********************
esc.ord = c('H1','H9','Sheff','T3ES','UVB01','SA01','oocyte')
i.esc.ord = unlist(lapply(esc.ord,grep,pt.name))

tis.e.mat = exprs.v[,c(i.lng.p,i.atrt.p,i.gbm.p,i.br.p,i.or.p,i.cn.p,i.lng.csc,i.atrt.csc,i.gbm.csc,i.br.csc,i.or.csc,i.cn.csc)]
esc.e.mat = exprs.v[,i.esc.ord]

t.out3 = cbind(t.out2,tis.e.mat,esc.e.mat)
t.out4 = cbind(tis.e.mat,esc.e.mat)
write.table(t.out3,"lvcom_mgX_data.txt",sep="\t", col.names =T, row.names=F, quote=F)




# ***************************************************************************************************************************************** #
# extract g.rxn 
# 
# **************************************************************************************************************************************** #

# ** extrac g.rxn ---- interaction network ------ from gene.Rxn ------ ** #
g.rxn <-  read.delim("hs_rxn_g.txt", sep="\t",header=T)
cat(c("The interaction data contains ", nrow(g.rxn), " rows, and ", ncol(g.rxn), "columns\n"))
col.grxn = c("id.a","id.b","pmid","txt","db")
colnames(g.rxn) = col.grxn

casi.r2 = ext.list.M2(nr.gid,g.rxn)		# *  rxns (nrow)  --> both of the interaction node are in the gene list
dim(casi.r2)

# -------------- remove self-directing links -------------- #							
gene.a =as.character(casi.r2[,1][drop=T])
sym.a = get.sym(gene.a,nr.gid,gid.sym)
gene.b =as.character(casi.r2[,2][drop=T])
sym.b = get.sym(gene.b,nr.gid,gid.sym)
casi.rxn = cbind(sym.a,sym.b,casi.r2)
ind = rep(T,length(sym.a))
ind[sym.a==sym.b]=F
casi.rxn2 = casi.rxn[ind,]			# 115


# *************************** remove Repetitive links *************************** #
tst.g = sprintf("%s_%s",casi.rxn2[,1],casi.rxn2[,2])
nr.grxn = removeRepeat(tst.g)
g.db = as.character(casi.rxn2[,7])
r.g.txt = sprintf("%s :: %s",casi.rxn2[,5],casi.rxn2[,6])
nw.out = c('ed-type','a.sym','b.sym','num-db','nw-db','txt')

for (i in 1:length(nr.grxn)) {
	pat = sprintf("^%s\\b$", nr.grxn[i])
	
	
	m = regexpr(pat,tst.g)
	p.db = paste(g.db[m>0],collapse="-")
	p.n.db = length(m[m>0])	
	p.txt = paste(r.g.txt[m>0],collapse=" ")
	
	p = unlist(strsplit(nr.grxn[i],"_"))
	p.a = p[1]
	p.b = p[2]
	p.ad = c("lvcom_grxn",p.a,p.b,p.n.db,p.db,p.txt)
	nw.out = rbind(nw.out,p.ad)
}
	
nw.out = nw.out[2:nrow(nw.out),]
colnames(nw.out) = c('ed-type','a.sym','b.sym','num-db','nw-db','txt')

a.gid = nr.gid[match(nw.out[,2],gid.sym)]
b.gid = nr.gid[match(nw.out[,3],gid.sym)]
nnw.out = cbind(nw.out[,1:3],a.gid,b.gid,nw.out[,4:6])

nw.grxn = nnw.out


# ***************************************************************************************************************************************** #
# merge IPA with g.rxn
# establish edge attribute
# start(link)target | class
# ***************************************************************************************************************************************** #

# ** ---------- read-in IPA network ----------- ** #
nw.ipa = read.delim('lvcom_icord_nw.txt',sep='\t',header=T,strip.white=T)
nw.ipa = cbind(nw.ipa,rep("IPA",nrow(nw.ipa)))

# ---------- ck and merge (1st) ------------- #
tst.g = sprintf("%s_%s",nw.grxn[,2],nw.grxn[,3])
tst.p = sprintf("%s_%s",nw.ipa[,2],nw.ipa[,3])
tst = c(tst.g,tst.p)
g.lg = nrow(nw.grxn)
nr.nw = removeRepeat(c(tst.g,tst.p))
ov.nw = outRepeat(c(tst.g,tst.p))
# merge
ed.tp = rep("lvcom_grxn-ipa",length(ov.nw))
ov.g.id = nw.grxn[match(ov.nw,tst.g),2:5]
ov.g.n.db = as.numeric(nw.grxn[match(ov.nw,tst.g),6])
ov.g.db = as.character(nw.grxn[match(ov.nw,tst.g),7])
ov.g.txt = as.character(nw.grxn[match(ov.nw,tst.g),8])

ov.p.n.db = as.numeric(nw.ipa[match(ov.nw,tst.p),6])
ov.p.db = as.character(nw.ipa[match(ov.nw,tst.p),7])
ov = cbind(ed.tp, ov.g.id, ov.g.n.db+ov.p.n.db, sprintf("%s-IPA.n%s",ov.g.db,ov.p.db),ov.g.txt)
write.table(ov,"lvcom_mg_Net.txt",sep="\t", col.names =F, row.names=F, quote=F)

# ---------- ck and merge (2nd) ------------- #
ind.p = match(tst.p,ov.nw,nomatch=0)
nw.nr.p = nw.ipa[ind.p==0,]
tt.p = sprintf("%s_%s",as.character(nw.nr.p[,3][drop=T]),as.character(nw.nr.p[,2][drop=T]))
ind.g = match(tst.g,ov.nw,nomatch=0)
nw.nr.g = nw.grxn[ind.g==0,]
tt.g = sprintf("%s_%s",as.character(nw.nr.g[,2][drop=T]),as.character(nw.nr.g[,3][drop=T]))
ov.tt = outRepeat(c(tt.p,tt.g))

# merge
ed.tp = rep("lvcom_grxn-ipa",length(ov.tt))
ov.g.id = nw.nr.g[match(ov.tt,tt.g),2:5]
ov.g.n.db = as.numeric(nw.nr.g[match(ov.tt,tt.g),6])
ov.g.db = as.character(nw.nr.g[match(ov.tt,tt.g),7])
ov.g.txt = as.character(nw.nr.g[match(ov.tt,tt.g),8])

ov.p.n.db = as.numeric(nw.nr.p[match(ov.tt,tt.p),6])
ov.p.db = as.character(nw.nr.p[match(ov.tt,tt.p),7])
ov.2 = cbind(ed.tp, ov.g.id, ov.g.n.db+ov.p.n.db, sprintf("%s-IPA.n%s",ov.g.db,ov.p.db),ov.g.txt)
write.table(ov.2,"lvcom_mg_Net.txt",sep="\t",append = T, col.names =F, row.names=F, quote=F)

# ---------- ck and write (3rd) ------------- #
ind.p = match(tt.p,ov.tt,nomatch=0)
nw.nr.p = nw.nr.p[ind.p==0,]
write.table(nw.nr.p,"lvcom_mg_Net.txt",sep="\t",append = T, col.names =F, row.names=F, quote=F)

ind.g = match(tt.g,ov.tt,nomatch=0)
nw.nr.g = nw.nr.g[ind.g==0,]
write.table(nw.nr.g,"lvcom_mg_Net.txt",sep="\t",append = T, col.names =F, row.names=F, quote=F)


# ***************************************************************************************************************************************** #
# extract exprs
# calculate corr
# ***************************************************************************************************************************************** #

# ** ---------- read-in final network ----------- ** #
nw.f = read.delim('lvcom_mg_Net.txt',sep='\t',header=F,strip.white=T)
a.id = as.character(nw.f[,4][drop=T])
a.pbs = g2pb.pbs[match(a.id,g2pb.gid)]
b.id = as.character(nw.f[,5][drop=T])
b.pbs = g2pb.pbs[match(b.id,g2pb.gid)]

nw.cor = c('match_probes','corr')

for(x in 1:nrow(nw.f)) {
	i.a = unlist(strsplit(a.pbs[x],","))
	i.b = unlist(strsplit(b.pbs[x],","))
	lg.a = length(i.a)
	lg.b = length(i.b)
	cat("nw.f of row:\t",x,"..........\n")
	
	ans.cor = c()
	ans.pb = c()
	for (i in 1:lg.a) {
			pba = i.a[i]
			x.a = exprs[match(pba,probe),ind.csc]
			for (j in 1:lg.b) {
				cat("this is :\t",i,"   ",j,"\n")
				pbb = i.b[j]
				x.b = exprs[match(pbb,probe),ind.csc]
				x.cor = cor(x.a,x.b)
				ans.cor = append(ans.cor,x.cor)
				ans.pb = append(ans.pb,sprintf("%s\t%s",pba,pbb))
			}
	}
	max.cor = ans.cor[ans.cor==max(ans.cor)]
	max.pbs = ans.pb[ans.cor==max(ans.cor)]
	nw.cor = rbind(nw.cor,c(max.pbs,max.cor))
}
nw.cor = nw.cor[2:nrow(nw.cor),]
ck = as.numeric(nw.cor[,2][drop=T])
ind.ck = rep(99,length(ck))
ind.ck[abs(ck)>=0.4] = 1
ind.ck[abs(ck)<0.4] = 0
length(ind.ck[ind.ck==0])
length(ind.ck[ind.ck==1])

nw.ff = cbind(nw.f,nw.cor,ind.ck)
write.table(nw.ff,"lvcom_mg_Net_cr40.txt",sep="\t",col.names =F, row.names=F, quote=F)		# complete network with corr
	
# ====================================================================================================================================
nw = nw.ff[ind.ck==1,]
# --------------- topology analysis ----------------- #
node.a = as.character(nw[,2][drop=T])
node.b = as.character(nw[,3][drop=T])
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

focality.ind = match(foc.g,label)   		
node.type[focality.ind] = 1
flow.ind = match(flow.g,label)				
node.type[flow.ind] = 2			
peri.ind = match(peri.g,label)								
node.type[peri.ind] = 3
node.sum = cbind(label,node.type,out.nbhood,out.foc.cp)		


# ---********** extract annotation *************---
# * t.out3 (previously processed)
# --------------------------------------------------
t.gsym = as.character(t.out3[,6][drop=T])
n.ind = match(label,t.gsym,nomatch=0)
t.dat = t.out3[n.ind,]
t.gid = as.character(t.dat[,4][drop=T])

n1.ind = match(label,a1.gsym,nomatch=0)					# only 12 genes mapped from original input
t.a1.dat = a1.dat[n1.ind,]
t.a1.gid = as.character(t.a1.dat[,10][drop=T])
t.x.dat = array("",dim=c(nrow(t.dat),8))
t.x.ind = match(t.a1.gid,t.gid,nomatch=0)
for (t in 1:length(t.x.ind)) {
	t.r = t.x.ind[t] 
	t.fill = t.a1.dat[t,c(2:6,11,7:8)]
	t.x.dat[t.r,1] = as.character(unlist(t.fill[1]))
	t.x.dat[t.r,2] = as.numeric(unlist(t.fill[2]))
	t.x.dat[t.r,3:6] = as.character(unlist(t.fill[3:6])[drop=T])
	t.x.dat[t.r,7:8] = as.numeric(unlist(t.fill[7:8])[drop=T])
}

n2.ind = match(a2.gid,t.gid,nomatch=0)					# one gene is pull up from 'grxn'-only
t.a2.dat = a2.dat[n2.ind,]
t.a2.gid = as.character(t.a2.dat[,7])
t.x2.dat = array("",dim=c(nrow(t.dat),7))
t.x2.ind = match(t.a2.gid,t.gid,nomatch=0)
for (t in 1:length(t.x2.ind)) {
	t.r = t.x2.ind[t] 
	t.fill = t.a2.dat[t,c(1:5,7:8)]
	t.x2.dat[t.r,1:2] = as.character(unlist(t.fill[1:2]))
	t.x2.dat[t.r,3] = as.numeric(unlist(t.fill[3]))
	t.x2.dat[t.r,4:5] = as.character(unlist(t.fill[4:5])[drop=T])
	t.x2.dat[t.r,6:7] = as.numeric(unlist(t.fill[6:7]))
}

node.out = cbind(node.sum,t.x.dat,t.x2.dat,t.dat[,1:9])

write.table(node.out,"lvcom_node_cr40.txt",sep="\t", col.names =T, row.names=F, quote=F)


