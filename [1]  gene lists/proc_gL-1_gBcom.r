# ----------------------------------------------------------------------------------------------------------------------------------------- #
# proc1 generate gene list gB_com
# ----------------------------------------------------------------------------------------------------------------------------------------- #

library('affyQCReport')
library('limma')
library('Mfuzz')
library('fdrtool')
source("sub_func.r")


# ----- annotation ----- # 										# --- annotation --- #
x = read.delim("hgU133p2_ann.txt", sep="\t",header=T)
probe = as.character(x[,1][drop = T])

ann.dat <-  read.delim("ca_probe.txt", sep="\t",header=F)
cat(c("The annotation data contains ", nrow(ann.dat), " rows, and ", ncol(ann.dat), "columns\n"))

ann.probe = trimWhiteSpace(as.character(ann.dat[,1][drop = T]))



# ------------------------- read-in pre-processed files --------------------- # 			
# quntSTD_ALL_z.txt
# -------------------------------------------------------------------- #			
data = read.table("quntSTD_ALL_z.txt",sep="\t",header=T,strip.white=T)
cat("reading file: \n")
cat("check data:	",dim(data),"\n")

pt = colnames(data)
pt.name =  gsub(".CEL", "", pt)
pt.name = gsub(".cel", "", pt.name)
pt.name = gsub("\\.", "_", pt.name)

# [1] "A549_cd133_1"   "A549_cd133_2"   "A549_cd133_3"   "A549_ir"        "A549_n1"        "A549_n2"        "A549_n3"        "A549_snail"     "ATRT_n"        
# [10] "ATRT_sf"        "BT_n"           "BT12_n"         "BT16_n"         "FaDu_bmi1"      "FaDu_p1"        "FaDu_p2"        "FaDu_twist"     "H1_gfp_1"      
# [19] "H1_o20_1"       "H1_o20_2"       "H1_o4_2"        "H9_1"           "H9_o20"         "H9_o4"          "H9_pool_1"      "H9_pool_2"      "H9_pool_3"     
# [28] "HT29_n1"        "HT29_n2"        "HT29_n3"        "HT29_n4"        "HT29_n5"        "HT29_sf1"       "HT29_sf2"       "HTB186_n"       "HTB186_sf"     
# [37] "MCF7_n1"        "MCF7_n2"        "MCF7_n3"        "MCF7_snail_6SA" "MCF7_snail"     "oocyte_1"       "oocyte_2"       "oocyte_3"       "PT1_ir1"       
# [46] "PT1_ir2"        "PT1_ir3"        "PT1_ir4"        "PT1_s1"         "PT1_s2"         "PT1_s3"         "SA01_1"         "SA01_2"         "SA01_3"        
# [55] "SAS_Luc"        "SAS_p_L"        "SAS_p_y"        "SAS_sf2_L"      "SAS_sf2_y"      "SAS_sf3_L"      "SAS_sf3_y"      "SAS_sf5_L"      "SAS_sf5_y"     
# [64] "Sheff4_1"       "Sheff4_2"       "Sheff4_3"       "SW480_n1"       "SW480_n2"       "SW480_snail_1"  "SW480_snail_2"  "T3ES_1"         "T3ES_2"        
# [73] "T3ES_3"         "U87mg_n"        "U87mg_sf"       "UVB01_1"        "UVB01_2"        "UVB01_3"   

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

# --- all genes in the g.List --- #
index = match(ann.probe,probe)			 	#13048
exprs.v = data.matrix(exprs[index,])
dim(exprs.v)


# --- ck exprs variation --- #
sd.p = apply(exprs.v[,ind.p],1,sd)
mean.p = apply(exprs.v[,ind.p],1,mean)
tq.p = sd.p/abs(mean.p)
ck.v.p = summary(tq.p)
ck.v.p
cut.p = sort(tq.p)[500]
cut.p
#cut.p = 0.5*(ck.v.p[1]+ck.v.p[2])
v.pb.p = ann.probe[tq.p<=cut.p]				# 
length(v.pb.p)

sd.csc = apply(exprs.v[,ind.csc],1,sd)
mean.csc = apply(exprs.v[,ind.csc],1,mean)
tq.csc = sd.csc/abs(mean.csc)
ck.v.csc = summary(tq.csc)
ck.v.csc
cut.csc = sort(tq.csc)[500]
cut.csc
#cut.csc = 0.5*(ck.v.csc[1]+ck.v.csc[2])
v.pb.csc = ann.probe[tq.csc<=cut.csc]		# 
length(v.pb.csc)

sd.Esc = apply(exprs.v[,ind.Esc],1,sd)
mean.Esc = apply(exprs.v[,ind.Esc],1,mean)
tq.Esc = sd.Esc/abs(mean.Esc)
ck.v.Esc = summary(tq.Esc)
ck.v.Esc
cut.Esc= sort(tq.Esc)[500]
cut.Esc
#cut.Esc = 0.5*(ck.v.Esc[1]+ck.v.Esc[2])
v.pb.Esc = ann.probe[tq.Esc<=cut.Esc]		# 
length(v.pb.Esc)


# ========================================================================================================================================= #
# * analysis / comparison

# * 1 * ---- ck overlapping area ---------------------------------------- *
w.csc.Esc = outRepeat(c(v.pb.csc,v.pb.Esc))			
length(w.csc.Esc)													# same (csc+Esc) 	= 77
w.p.csc = outRepeat(c(v.pb.p,v.pb.csc))		
length(w.p.csc)														# same (csc+p) 		= 196
w.p.Esc = outRepeat(c(v.pb.p,v.pb.Esc))		
length(w.p.Esc)														# same (p+Esc)		= 91
w.ALL3 = outRepeat(c(w.p.csc,w.csc.Esc)	)	
length(w.ALL3)														# same (p+csc+Esc) 	= 49


# * 2 * ---- overlap specific area (2) ---------------------------------- *
i.v1 = rep(T,length(w.p.csc))
i.v1[match(w.ALL3,w.p.csc)]=F
lv_pCsc_nEsc = w.p.csc[i.v1]			
length(lv_pCsc_nEsc)												# 147

i.v2 = rep(T,length(w.p.Esc))
i.v2[match(w.ALL3,w.p.Esc)]=F
lv_pEsc_nCsc = w.p.Esc[i.v2]			
length(lv_pEsc_nCsc)												# 42 

i.v3 = rep(T,length(w.csc.Esc))
i.v3[match(w.ALL3,w.csc.Esc)]=F
lv_cscEsc_nP = w.csc.Esc[i.v3]			
length(lv_cscEsc_nP)												# 28

# * 3 * ---- only specific area (1) ------------------------------------- *
i.v4 = rep(T,length(v.pb.Esc))
i.v4[match(w.p.Esc,v.pb.Esc,nomatch=0)]=F
i.v4[match(w.csc.Esc,v.pb.Esc,nomatch=0)]=F
lv_Esc_nPnCsc = v.pb.Esc[i.v4]							
length(lv_Esc_nPnCsc)												# 381   

i.v5 = rep(T,length(v.pb.p))
i.v5[match(w.p.Esc,v.pb.p,nomatch=0)]=F
i.v5[match(w.p.csc,v.pb.p,nomatch=0)]=F
lv_p_nEsc_nCsc = v.pb.p[i.v5]							
length(lv_p_nEsc_nCsc)												# 262

i.v6 = rep(T,length(v.pb.csc))
i.v6[match(w.p.csc,v.pb.csc,nomatch=0)]=F
i.v6[match(w.csc.Esc,v.pb.csc,nomatch=0)]=F
lv_csc_nPnEsc = v.pb.csc[i.v6]							
length(lv_csc_nPnEsc)												# 276

# * 4 * ---- define specific g.lv.* ---------------------- ***
g.lv1.p = 			lv_p_nEsc_nCsc			# 262
g.lv2.Esc = 		lv_Esc_nPnCsc			# 381		
g.lv3.csc = 		lv_csc_nPnEsc			# 276
g.lv4.p.Esc	= 		lv_pEsc_nCsc			# 42
g.lv5.p.csc	=		lv_pCsc_nEsc			# 147
g.lv6.csc.Esc = 	lv_cscEsc_nP			# 28
g.lv7.ALL = 		w.ALL3					# 49






# ========================================================================================================================================= #
# * 
# -------------------- Diff exprs probes ----------------------- #
set.a <-  read.delim("gA_tCnd.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(set.a), " rows, and ", ncol(set.a), "columns\n"))

i.pb = as.character(set.a[,1][drop=T])
i.com = as.numeric(set.a[,5][drop=T])
i.com.pb = i.pb[i.com>=2]

# -- oral lung gbm atrt br cn snail sf
i.or = as.numeric(set.a[,6][drop=T])
i.or.pb = i.pb[i.or==1]
i.lg = as.numeric(set.a[,7][drop=T])
i.lg.pb = i.pb[i.lg==1]
i.gbm = as.numeric(set.a[,8][drop=T])
i.gbm.pb = i.pb[i.gbm==1]
i.atrt = as.numeric(set.a[,9][drop=T])
i.atrt.pb = i.pb[i.atrt==1]
i.br = as.numeric(set.a[,10][drop=T])
i.br.pb = i.pb[i.br==1]
i.cn = as.numeric(set.a[,11][drop=T])
i.cn.pb = i.pb[i.cn==1]
i.snail = as.numeric(set.a[,12][drop=T])
i.snail.pb = i.pb[i.snail==1]
i.sf = as.numeric(set.a[,13][drop=T])
i.sf.pb = i.pb[i.sf==1]


# --- overlap --- #
set.a.pb = i.pb
ov.1 = ck2.same(set.a.pb,g.lv1.p)
length(ov.1)										# 25
ov.2 = ck2.same(set.a.pb,g.lv2.Esc)
length(ov.2)										# 21
ov.3 = ck2.same(set.a.pb,g.lv3.csc)
length(ov.3)										# 25
ov.4 = ck2.same(set.a.pb,g.lv4.p.Esc)
length(ov.4)										# 4
ov.5 = ck2.same(set.a.pb,g.lv5.p.csc)
length(ov.5)										# 10
ov.6 = ck2.same(set.a.pb,g.lv6.csc.Esc)
length(ov.6)										# 4
ov.7 = ck2.same(set.a.pb,g.lv7.ALL)
length(ov.7)										# 5

# * total overlap
i.v7 = rep(T,length(set.a.pb))						# set.a =  789 = 695 + 94
i.v7[match(v.pb.p,set.a.pb,nomatch=0)]=F
i.v7[match(v.pb.csc,set.a.pb,nomatch=0)]=F
i.v7[match(v.pb.Esc,set.a.pb,nomatch=0)]=F
nov_set.a = set.a.pb[i.v7]							
length(nov_set.a)									# 695
lv_set.a = c(ov.1,ov.2,ov.3,ov.4,ov.5,ov.6,ov.7)
length(lv_set.a)									# 94

# * lv_csc not overlap with set.a *
i.v8 = rep(T,length(v.pb.csc))
i.v8[match(lv_set.a,v.pb.csc,nomatch=0)]=F
lv_csc_nSetA = v.pb.csc[i.v8]							
length(lv_csc_nSetA)								# 456


# ********************************** index matrix ****************************************
# list: 
# v.pb.csc, set.a.pb
# lv_set.a, nov_set.a, lv_csc_nSetA
# g.lv7.ALL, g.lv6.csc.Esc, g.lv5.p.csc, g.lv3.csc
# tcnd (ind)
# ---------------------------------------------------------------------------------------- 

gnr.pb = removeRepeat(c(set.a.pb,v.pb.csc))
lg = length(gnr.pb)

# tissue-cond
tcnd = rep("x",lg)
tcnd[match(i.or.pb,gnr.pb)] = "or"
tcnd[match(i.lg.pb,gnr.pb)] = "lg"
tcnd[match(i.gbm.pb,gnr.pb)] = "gbm"
tcnd[match(i.atrt.pb,gnr.pb)] = "atrt"
tcnd[match(i.br.pb,gnr.pb)] = "br"
tcnd[match(i.cn.pb,gnr.pb)] = "cn"
tcnd[match(i.snail.pb,gnr.pb)] = "snail"
tcnd[match(i.sf.pb,gnr.pb)] = "sf"
tcnd[match(i.com.pb,gnr.pb)] = "com"

# ======================== corr ============================ #
# already calculated -- not run, very time-consuming
# ---------------------------------------------------------- #
i.cor = match(gnr.pb,probe)			 	
exprs.v = data.matrix(exprs[i.cor,])
dim(exprs.v)	# 1245 78

# --- Csc --- #
lg.pb = length(gnr.pb)
cor.mat.csc = array(NA,dim = c(lg.pb,lg.pb))

for (i in 1:lg.pb) {
	for (j in 1:lg.pb) {
		exprs.g.a = exprs.v[i,ind.csc]
		exprs.g.b = exprs.v[j,ind.csc]
		cor.mat.csc[i,j] = cor(exprs.g.a,exprs.g.b)
	}
}
cor.mean.csc = apply(abs(cor.mat.csc),1,mean)

# --- P --- #
cor.mat.p = array(NA,dim = c(lg.pb,lg.pb))

for (i in 1:lg.pb) {
	for (j in 1:lg.pb) {
		exprs.g.a = exprs.v[i,ind.p]
		exprs.g.b = exprs.v[j,ind.p]
		cor.mat.p[i,j] = cor(exprs.g.a,exprs.g.b)
	}
}
cor.mean.p = apply(abs(cor.mat.p),1,mean)

# --- Esc --- #
cor.mat.Esc = array(NA,dim = c(lg.pb,lg.pb))

for (i in 1:lg.pb) {
	for (j in 1:lg.pb) {
		exprs.g.a = exprs.v[i,ind.Esc]
		exprs.g.b = exprs.v[j,ind.Esc]
		cor.mat.Esc[i,j] = cor(exprs.g.a,exprs.g.b)
	}
}
cor.mean.Esc = apply(abs(cor.mat.Esc),1,mean)

summary(cor.mean.csc)
summary(cor.mean.p)
summary(cor.mean.Esc)
cor.out = cbind(cor.mean.csc,cor.mean.p,cor.mean.Esc)
write.table(cor.out,"g_x_cor.txt",sep="\t",quote=F,col.names=T,row.names=F)

# ============================================================ #


# ******************** #
# *** index matrix *** #
# ******************** #
ind.m = matrix(0,nrow=lg,ncol=10)
ind.m[match(v.pb.csc,gnr.pb),1]=1
ind.m[match(set.a.pb,gnr.pb),2]=1
ind.m[match(lv_set.a,gnr.pb),3]=1
ind.m[match(nov_set.a,gnr.pb),4]=1
ind.m[match(lv_csc_nSetA,gnr.pb),5]=1
ind.m[match(g.lv7.ALL,gnr.pb),6]=1
ind.m[match(g.lv6.csc.Esc,gnr.pb),7]=1
ind.m[match(g.lv5.p.csc,gnr.pb),8]=1
ind.m[match(g.lv3.csc,gnr.pb),9]=1
ind.m[,10]= tcnd

colnames(ind.m) = c("v.pb.csc","set.a.pb","lv_set.a","nov_set.a","lv_csc_nSetA",
                    "g.lv7.ALL","g.lv6.csc.Esc","g.lv5.p.csc","g.lv3.csc","tcnd")

# ----- annotation ----- # 										# --- annotation --- #
ann.dat <-  read.delim("hgu133p2_ann.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(ann.dat), " rows, and ", ncol(ann.dat), "columns\n"))
a.pb = as.character(ann.dat[,1][drop = T])

ind.ann = ann.dat[match(gnr.pb,a.pb),]

out = cbind(ind.ann,ind.m)
write.table(out,"g_Com.txt",sep="\t",quote=F,col.names=T,row.names=F)




