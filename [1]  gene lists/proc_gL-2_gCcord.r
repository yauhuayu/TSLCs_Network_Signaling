# ----------------------------------------------------------------------------------------------------------------------------------------- #
# proc2 generate gene list gC_final_cord
# ----------------------------------------------------------------------------------------------------------------------------------------- #

library('affyQCReport')
library('limma')
library('Mfuzz')
library('fdrtool')
source("sub_func.r")


# ----- annotation ----- # 																							# --- annotation --- #
ann.dat <-  read.delim("gB_Com.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(ann.dat), " rows, and ", ncol(ann.dat), "columns\n"))
# ann.probe	fdf.pCsc	fdf.EscCsc	g.id.num	out.gid	g.sym	g.chr	g.omim	n.ann.tsp	g.name	
# [11] v.pb.csc	set.a.pb	lv_set.a	nov_set.a	lv_csc_nSetA	g.lv7.ALL	g.lv6.csc.Esc	g.lv5.p.csc	g.lv3.csc	tcnd	
# [21] cor.mean.csc	cor.mean.p	cor.mean.Esc
ann.probe = as.character(ann.dat[,1][drop = T])
g.name = as.character(ann.dat[,10][drop = T])
g.sym = as.character(ann.dat[,6][drop = T])
g.id = as.character(ann.dat[,5][drop = T])


# ------------------------- read-in pre-processed files --------------------------------------------------------------- read-in ---------- # 			
# quntSTD_ALL.txt
# ---------------------------------------------------------------------------------------------------------------------------------------- #			
data = read.table("quntSTD_ALL_z.txt",sep="\t",header=T,strip.white=T)
cat("reading file: \n")
cat("check data:	",dim(data),"\n")
pt = colnames(data)
pt.name =  gsub(".CEL", "", pt)
pt.name = gsub(".cel", "", pt.name)
pt.name = gsub("\\.", "_", pt.name)
probe = rownames(data)

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

# --- all genes in the g.Com --- #
index = match(ann.probe,probe)			 	#
exprs.v = data.matrix(exprs[index,])

# -------- ck mean fdf ----------- missing in proc_1 file ---------- #
ck.p = exprs.v[,ind.p]
ck.csc = exprs.v[,ind.csc]
ck.Esc = exprs.v[,ind.Esc]

ck.mean.p = apply(ck.p,1,mean)
ck.mean.csc = apply(ck.csc,1,mean)
ck.mean.Esc = apply(ck.Esc,1,mean)

ck.fdf.pCsc = ck.mean.csc - ck.mean.p
ck.fdf.EscCsc = ck.mean.Esc - ck.mean.csc
ck.out = cbind(ann.probe,ck.fdf.pCsc,ck.fdf.EscCsc)

write.table(ck.out,"ck_X_data.txt",sep="\t", col.names =T, row.names=F, quote=F)



# ------------------------- tissue type marker -------------------------------------------------------------------------------- # 			
#  col: 
#  1-4: 	v.pb.csc 	set.a.pb  	tcnd 		lv_x1_a 
#  5-9:		lv_x2_c 	glv_1_all 	glv_2_ec 	glv_3_pc 	glv_4_c 
#  10-13:	glv_5_a     cor_c     	cor_p     	cor_e
# ----------------------------------------------------------------------------------------------------------------------------- #
ind.m = ann.dat[,11:ncol(ann.dat)]  # 7:19

# [1] "A549_cd133_1"   "A549_cd133_2"   "A549_cd133_3"   "A549_ir"        "A549_n1"        "A549_n2"        "A549_n3"        "A549_snail"     "ATRT_n"        
# [10] "ATRT_sf"        "BT_n"           "BT12_n"         "BT16_n"         "FaDu_bmi1"      "FaDu_p1"        "FaDu_p2"        "FaDu_twist"     "H1_gfp_1"      
# [19] "H1_o20_1"       "H1_o20_2"       "H1_o4_2"        "H9_1"           "H9_o20"         "H9_o4"          "H9_pool_1"      "H9_pool_2"      "H9_pool_3"     
# [28] "HT29_n1"        "HT29_n2"        "HT29_n3"        "HT29_n4"        "HT29_n5"        "HT29_sf1"       "HT29_sf2"       "HTB186_n"       "HTB186_sf"     
# [37] "MCF7_n1"        "MCF7_n2"        "MCF7_n3"        "MCF7_snail_6SA" "MCF7_snail"     "oocyte_1"       "oocyte_2"       "oocyte_3"       "PT1_ir1"       
# [46] "PT1_ir2"        "PT1_ir3"        "PT1_ir4"        "PT1_s1"         "PT1_s2"         "PT1_s3"         "SA01_1"         "SA01_2"         "SA01_3"        
# [55] "SAS_Luc"        "SAS_p_L"        "SAS_p_y"        "SAS_sf2_L"      "SAS_sf2_y"      "SAS_sf3_L"      "SAS_sf3_y"      "SAS_sf5_L"      "SAS_sf5_y"     
# [64] "Sheff4_1"       "Sheff4_2"       "Sheff4_3"       "SW480_n1"       "SW480_n2"       "SW480_snail_1"  "SW480_snail_2"  "T3ES_1"         "T3ES_2"        
# [73] "T3ES_3"         "U87mg_n"        "U87mg_sf"       "UVB01_1"        "UVB01_2"        "UVB01_3"   

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


# ------------------------------------------------------------------------------------------------ #
# check direction
# ------------------------------------------------------------------------------------------------ #
tis.cond = colnames(ind.m)
for(h in c(3,4,6:9)){
	y = h
	i.y = ind.m[,y]
	i.p = i.sf.p									# ** need to change tis-type per loop ** #
	i.c = i.sf.csc

	lg.p = length(i.p)
	lg.c = length(i.c)
	lg.pbs = length(ann.probe[i.y==1])
	epx.p = exprs.v[i.y==1,i.p]
	epx.c = exprs.v[i.y==1,i.c]

	out = c("total","up","dn","major")

	for (nr in 1:lg.pbs) {
		#nr = 1
		ans = c()
		for (i in 1:lg.p) {
				x.p = epx.p[nr,i]
				for (j in 1:lg.c) {
					#cat("this is :\t",i,"   ",j,"\n")
					x.c = epx.c[nr,j]
					x.df = x.c-x.p
					ans = append(ans,x.df)
				}
		}
		tot = length(ans)
		up = length(ans[ans>0])
		dn = length(ans[ans<0])
		majr = max(up,dn)/tot
		out = rbind(out,c(tot,up,dn,majr))
	}
		
	out = data.matrix(out[2:nrow(out),])
	ck.maj = as.numeric(out[,4])
	out = cbind(ann.dat[i.y==1,c(1:7,9,10,20:23)],out)
	out2 = out[ck.maj==1,]
	dim(out2)
	out2 = cbind(rep(tis.cond[y],nrow(out2)),out2)
	write.table(out2,"out_test.txt",sep="\t", append = F,row.names=F,col.names=F,quote=F)

}   	# type-cond

tst.pb <-  read.delim("out_test.txt", sep="\t",header=F)
cat(c("The annotation data contains ", nrow(tst.pb), " rows, and ", ncol(tst.pb), "columns\n"))

catg = as.character(tst.pb[,1][drop = T])
tis.pb = as.character(tst.pb[,2][drop = T])
xnr.pb = removeRepeat(tis.pb)
xnr.catg = rep(NA,length(xnr.pb))
xnr.lg = rep(0,length(xnr.pb))
for(w in 1:length(xnr.pb)) {
	xnr.catg[w] = paste(catg[grep(xnr.pb[w],tis.pb)],collapse = " ")
	xnr.lg[w] = length(catg[grep(xnr.pb[w],tis.pb)])
}
f.out = cbind(xnr.pb,xnr.lg,xnr.catg,rep("sf",length(xnr.pb)))				 			# ** change tis-type ** #

write.table(f.out,"out_concord_pbs.txt",sep="\t", append = T,row.names=F,col.names=F,quote=F)
	
# ===================================================================================== #
# ck final 
# ===================================================================================== #
# change/merge final 'out_concord_pbs.txt' to 'out_tis_pbs.txt' 
# after looping all comparisons -- lung, atrt, oral, colon, breast, colon, snail, sf --
tst.pb <-  read.delim("out_tis_pbs.txt", sep="\t",header=F)					
cat(c("The annotation data contains ", nrow(tst.pb), " rows, and ", ncol(tst.pb), "columns\n"))

tis.pb = as.character(tst.pb[,1][drop = T])
tis.cnt = as.numeric(tst.pb[,2][drop = T])
tis.sdb = as.character(tst.pb[,3][drop = T])
tis = as.character(tst.pb[,4][drop = T])
xnr.pb = removeRepeat(tis.pb)
xnr.lg = rep(0,length(xnr.pb))
xnr.tis = rep(NA,length(xnr.pb))
xnr.tis.cnt = rep(NA,length(xnr.pb))
xnr.tis.sdb = rep(NA,length(xnr.pb))

for(w in 1:length(xnr.pb)) {
	xnr.lg[w] = length(tis[grep(xnr.pb[w],tis.pb)])
	xnr.tis[w]= paste(tis[grep(xnr.pb[w],tis.pb)],collapse = " ")
	xnr.tis.cnt[w]= paste(tis.cnt[grep(xnr.pb[w],tis.pb)],collapse = ",")
	xnr.tis.sdb[w]= paste(tis.sdb[grep(xnr.pb[w],tis.pb)],collapse = ",")
}
f.out = cbind(xnr.pb,xnr.lg,xnr.tis,xnr.tis.cnt,xnr.tis.sdb,ann.dat[match(xnr.pb,ann.probe),c(1:7,9,10,20:23)])	 		
write.table(f.out,"out_cord_FINAL.txt",sep="\t", append = T,row.names=F,col.names=F,quote=F)
	