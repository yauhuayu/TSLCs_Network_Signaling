# ----------------------------------------------------------------------------------------------------------------------------------------------- #
# proc5_lvcom_state.r
# network state analysis of (p :: csc :: Esc)
# ----------------------------------------------------------------------------------------------------------------------------------------------- #

library('limma')
library('Matching')																										    # ks.boot
source("sub_func.r")
source("method.r")


# ** ------------- (1) process core_network (summary) --------------------------------------------------------------------- *** (1) *** --- #
# * read in exprs data 
# ** -------------------------------------------------------------------------------------------------------------------------------------- #
n.dat <-  read.delim("lvcom_node_cr40.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(n.dat), " rows, and ", ncol(n.dat), "columns\n"))
n.gsym= as.character(n.dat[,1][drop = T])
n.type = as.numeric(n.dat[,2][drop=T])
n.nb = as.numeric(n.dat[,3][drop=T])
n.foc = as.numeric(n.dat[,4][drop=T])
n.pb = as.character(n.dat[,7][drop = T])
n.gname = as.character(n.dat[,14][drop = T])


# ** ------------- (2) consensus network signal dynamcis analysis -----------------------------------------------       ---- *** (2a) *** ---- 
# read - in exprs
# ----------------
exprs = read.table("lvcom_mg_data_z2.txt", sep="\t",header=T)
e.sym = as.character(exprs[,6][drop=T])
i.net = match(n.gsym, e.sym)
cbind(n.gsym,exprs[i.net,c(1:8)])[1:5,]			  # for ck 

epx.net = exprs[i.net,10:ncol(exprs)]
sample = colnames(epx.net)
i.p = 1:27
i.csc = 28:(27+26)				# 53
i.esc = 54:78

# [1] "A549_n1"        "A549_n2"        "A549_n3"        "ATRT_n"         "BT_n"           "BT12_n"         "BT16_n"         "HTB186_n"       "PT1_s1"        
# [10] "PT1_s2"         "PT1_s3"         "U87mg_n"        "MCF7_n1"        "MCF7_n2"        "MCF7_n3"        "FaDu_p1"        "FaDu_p2"        "SAS_Luc"       
# [19] "SAS_p_L"        "SAS_p_y"        "HT29_n1"        "HT29_n2"        "HT29_n3"        "HT29_n4"        "HT29_n5"        "SW480_n1"       "SW480_n2"      
# [28] "A549_cd133_1"   "A549_cd133_2"   "A549_cd133_3"   "A549_ir"        "A549_snail"     "ATRT_sf"        "HTB186_sf"      "PT1_ir1"        "PT1_ir2"       
# [37] "PT1_ir3"        "PT1_ir4"        "U87mg_sf"       "MCF7_snail_6SA" "MCF7_snail"     "FaDu_bmi1"      "FaDu_twist"     "SAS_sf2_L"      "SAS_sf2_y"     
# [46] "SAS_sf3_L"      "SAS_sf3_y"      "SAS_sf5_L"      "SAS_sf5_y"      "HT29_sf1"       "HT29_sf2"       "SW480_snail_1"  "SW480_snail_2"  "H1_gfp_1"      
# [55] "H1_o20_1"       "H1_o20_2"       "H1_o4_2"        "H9_1"           "H9_o20"         "H9_o4"          "H9_pool_1"      "H9_pool_2"      "H9_pool_3"     
# [64] "Sheff4_1"       "Sheff4_2"       "Sheff4_3"       "T3ES_1"         "T3ES_2"         "T3ES_3"         "UVB01_1"        "UVB01_2"        "UVB01_3"       
# [73] "SA01_1"         "SA01_2"         "SA01_3"         "oocyte_1"       "oocyte_2"       "oocyte_3"      


# ===== ------------ ===== #
# grouping of node																										---- *** (2b) *** ----
# ===== ------------ ===== #
i.foc = grep(1,n.type)
i.flow = grep(2,n.type)
i.peri = grep(3,n.type)
cat("total number of genes:\t",length(n.gsym),"\n","inter - \t",length(i.foc),"\t intra - \t",length(i.flow),"\n")


# ----------------------------------------------------------------------------------- * re-ordering * ----------------  ---- *** (2c) *** ---- 
epx.ng = epx.net[c(i.foc,i.flow,i.peri),]																	

wt.0.nb =  c(rep(0,length(i.foc)),rep(1,length(i.flow)),rep(0,length(i.peri)))
wt.0.f =  c(rep(1,length(i.foc)),rep(0,length(i.flow)),rep(0,length(i.peri)))

wt.ab.nb = c(n.nb[c(i.foc,i.flow)],rep(0,length(i.peri)))
wt.ab.f = c(n.foc[i.foc],n.nb[i.flow],rep(0,length(i.peri)))

wt.a.nb = c(rep(0,length(i.foc)),n.nb[i.flow],rep(0,length(i.peri)))

wt.b.nb = c(n.nb[i.foc],rep(0,length(i.flow)),rep(0,length(i.peri)))
wt.b.f = c(n.foc[i.foc],rep(0,length(i.flow)),rep(0,length(i.peri)))

# ===================== calculate prediction  ============================ 						------------------------------ *** (3) *** ----
# ** epx **
pred.ab.epx = apply(epx.ng,2,sum)/nrow(epx.ng)
pred.a.epx = apply(epx.ng*wt.0.nb,2,sum)/nrow(epx.ng)
pred.b.epx = apply(epx.ng*wt.0.f,2,sum)/nrow(epx.ng)

# ** mag **
pred.ab.nb.mag.x = abs(epx.ng)*wt.ab.nb
pred.ab.nb.mag = apply(pred.ab.nb.mag.x,2,sum)/nrow(epx.ng)
pred.ab.f.mag.x = abs(epx.ng)*wt.ab.f
pred.ab.f.mag = apply(pred.ab.f.mag.x,2,sum)/nrow(epx.ng)

pred.a.nb.mag.x = abs(epx.ng)*wt.a.nb
pred.a.nb.mag = apply(pred.a.nb.mag.x,2,sum)/nrow(epx.ng)
pred.b.nb.mag.x = abs(epx.ng)*wt.b.nb
pred.b.nb.mag = apply(pred.b.nb.mag.x,2,sum)/nrow(epx.ng)
pred.b.f.mag.x = abs(epx.ng)*wt.b.f
pred.b.f.mag = apply(pred.b.f.mag.x,2,sum)/nrow(epx.ng)

# ** spec **
pred.ab.nb.spec.1 = apply(abs(epx.ng),2,k.b.wt.dis.G,wt.ab.nb)/(nrow(epx.ng)-1)
pred.ab.nb.spec.2 = apply(abs(epx.ng),2,k.b.wt.dis.G.contra,wt.ab.nb)/(nrow(epx.ng)-1)
pred.ab.nb.spec.x = 0.5*(pred.ab.nb.spec.1 + pred.ab.nb.spec.2)
pred.ab.nb.spec = apply(pred.ab.nb.spec.x,2,sum)/nrow(epx.ng)

pred.ab.f.spec.1 = apply(abs(epx.ng),2,k.b.wt.dis.G,wt.ab.f)/(nrow(epx.ng)-1)
pred.ab.f.spec.2 = apply(abs(epx.ng),2,k.b.wt.dis.G.contra,wt.ab.f)/(nrow(epx.ng)-1)
pred.ab.f.spec.x = 0.5*(pred.ab.f.spec.1 + pred.ab.f.spec.2)
pred.ab.f.spec = apply(pred.ab.f.spec.x,2,sum)/nrow(epx.ng)

pred.a.nb.spec.1 = apply(abs(epx.ng),2,k.b.wt.dis.G,wt.a.nb)/(nrow(epx.ng)-1)
pred.a.nb.spec.2 = apply(abs(epx.ng),2,k.b.wt.dis.G.contra,wt.a.nb)/(nrow(epx.ng)-1)
pred.a.nb.spec.x = 0.5*(pred.a.nb.spec.1 + pred.a.nb.spec.2)
pred.a.nb.spec = apply(pred.a.nb.spec.x,2,sum)/nrow(epx.ng)

pred.b.nb.spec.1 = apply(abs(epx.ng),2,k.b.wt.dis.G,wt.b.nb)/(nrow(epx.ng)-1)
pred.b.nb.spec.2 = apply(abs(epx.ng),2,k.b.wt.dis.G.contra,wt.b.nb)/(nrow(epx.ng)-1)
pred.b.nb.spec.x = 0.5*(pred.b.nb.spec.1 + pred.b.nb.spec.2)
pred.b.nb.spec = apply(pred.b.nb.spec.x,2,sum)/nrow(epx.ng)

pred.b.f.spec.1 = apply(abs(epx.ng),2,k.b.wt.dis.G,wt.b.f)/(nrow(epx.ng)-1)
pred.b.f.spec.2 = apply(abs(epx.ng),2,k.b.wt.dis.G.contra,wt.b.f)/(nrow(epx.ng)-1)
pred.b.f.spec.x = 0.5*(pred.b.f.spec.1 + pred.b.f.spec.2)
pred.b.f.spec = apply(pred.b.f.spec.x,2,sum)/nrow(epx.ng)

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


# ===============================================================================      # ** ------  plotting dynamics --- **
# * color indexing 
# lung - orange,brown
# atrt - lightgreen, green
# gbm - lightblue, blue
# br - pink, red
# or - grey, darkgrey
# cn - cyan, navy
# esc - yellow / oocyte - black

# * K-S test (bootstrap version)
q.type = c('P vs T','T vs E','P vs E')
q.out = q.type
w.all = list(pred.ab.epx,pred.a.epx,pred.b.epx,pred.ab.nb.mag,pred.ab.f.mag,pred.a.nb.mag,pred.b.nb.mag,pred.b.f.mag,pred.ab.nb.spec,pred.ab.f.spec,
pred.a.nb.spec,pred.b.nb.spec,pred.b.f.spec,pred.ab.nb.score,pred.ab.f.score,pred.a.nb.score,pred.b.nb.score,pred.b.f.score)

for (q in 1:length(w.all)) {
	w = w.all[[q]]
	w.p = w[i.p]
	w.csc = w[i.csc]
	w.esc = w[i.esc]
	ks.ck.pv = unlist(c(ks.boot(w.p,w.csc)$ks[2],ks.boot(w.csc,w.esc)$ks[2],ks.boot(w.p,w.esc)$ks[2]))
	cat("PV for KS-test:\n",ks.ck.pv[1],"(P vs T)\n",ks.ck.pv[2],"(T vs E)\n",ks.ck.pv[3],"(P vs E)\n\n")
	q.out = cbind(q.out,ks.ck.pv)
}

ck.q.col = c('ab.epx','a.epx','b.epx','ab.nb.mag','ab.f.mag','a.nb.mag','b.nb.mag','b.f.mag','ab.nb.spec','ab.f.spec',
'a.nb.spec','b.nb.spec','b.f.spec','ab.nb.SNR','ab.f.SNR','a.nb.SNR','b.nb.SNR','b.f.SNR')
rownames(q.out) = q.type
q.out = q.out[,2:ncol(q.out)]
colnames(q.out) = ck.q.col
# write.table(q.out,'test_ks_core.txt',sep="\t", col.names =T, row.names=F, quote=F)


w.col = c(rep('orange',3),rep('lightgreen',5),rep('lightblue',4),rep('pink',3),rep('grey',5),rep('cyan',7),rep('brown',5),rep('green',2),rep('blue',5),rep('red',2),rep('darkgrey',8),rep('navy',4),rep('yellow',22),rep('black',3))
w.sig = list(pred.b.epx,pred.b.nb.score)
w.title = c('Expression','SNR (wt=degree)')
q = 2
w = w.sig[[q]]
w.p = w[i.p]
w.csc = w[i.csc]
w.esc = w[i.esc]

plot(w.p,col = w.col[i.p],pch = 1, ylim = c(min(w)-1,max(w)+1.5),xlim =c(1,34), main=sprintf("Inter-modular hubs %s",w.title[q]),xlab="",ylab="dynamic range")
lines(w.p,col='green')
points(w.csc,col = w.col[i.csc],pch = 2)
lines(w.csc,col='red')
points(w.esc,col = w.col[i.esc],pch = 3)
lines(w.esc,col='orange')
t.ck.pv = c(t.test(w.p,w.csc)$p.value,t.test(w.csc,w.esc)$p.value,t.test(w.p,w.esc)$p.value)
cat("PV for T-test:\n",t.ck.pv[1],"(P vs CSC)\n",t.ck.pv[2],"(CSC vs ESC)\n",t.ck.pv[3],"(P vs ESC)\n\n")
ks.ck.pv = unlist(c(ks.boot(w.p,w.csc)$ks[2],ks.boot(w.csc,w.esc)$ks[2],ks.boot(w.p,w.esc)$ks[2]))
cat("PV for KS-test:\n",ks.ck.pv[1],"(P vs CSC)\n",ks.ck.pv[2],"(CSC vs ESC)\n",ks.ck.pv[3],"(P vs ESC)\n\n")

pvv = ks.ck.pv
txt.pCsc = sprintf("P vs T: pv = %.3f", pvv[1])
txt.cscEsc = sprintf("T vs E: pv = %.3f", pvv[2])
txt.pEsc = sprintf("P vs E: pv = %.3f", pvv[3])
txt = c(txt.pCsc,txt.cscEsc,txt.pEsc)
legend(0,max(w)+1.5,txt,cex=0.65,text.col = 'black', bg ='white',bty="n")
legend(1,max(w)+0.75,c("PTCs","TSLCs","hESCs"),cex=0.65,col=c('green','red','orange'),lty=1,text.col = 'black', bg ='white',bty="n")

legend(10,max(w)+1.5,sample,cex = 0.38,col = w.col,pch = c(rep(1,length(i.p)),rep(2,length(i.csc)),rep(3,length(i.esc))),ncol = 5,text.col='black',bg ='white',bty="n")



# ================= box plot ==================== #
#par(mfrow=c(1,3))
q = 2
w = w.sig[[q]]
w.p = w[i.p]
w.csc = w[i.csc]
w.esc = w[i.esc]
t.ck.pv = c(t.test(w.p,w.csc)$p.value,t.test(w.csc,w.esc)$p.value,t.test(w.p,w.esc)$p.value)

boxplot(w.p,w.csc,w.esc,ylim =c(min(w)-0.4,max(w)+0.4),col = c('green','red','orange'),main=w.title[q])
pvv = t.ck.pv
txt.pCsc = sprintf("PTCs vs TSLCs: pv = %.3f", pvv[1])
txt.cscEsc = sprintf("TSLCs vs hESCs: pv = %.3f", pvv[2])
txt.pEsc = sprintf("PTCs vs hESCs: pv = %.3f", pvv[3])
txt = c(txt.pCsc,txt.cscEsc,txt.pEsc)
legend(1,max(w)+0.4,txt,cex=0.75,text.col = 'black', bg ='white',bty="n")


# ============================================================================================ ** by gene ** ========================
n.label = n.gsym[c(i.foc,i.flow,i.peri)]
q.type = c('P vs T','T vs E','P vs E')
q.out = q.type

for (q in 1:length(n.label)) {
	#gind = match('DNMT3A',n.label)
	gind = q
	w =  epx.ng
	w.p = unlist(w[gind,i.p])
	w.csc = unlist(w[gind,i.csc])
	w.esc = unlist(w[gind,i.esc])
	w.ck = c(w.p,w.csc,w.esc)
	ks.ck.pv = unlist(c(ks.boot(w.p,w.csc)$ks[2],ks.boot(w.csc,w.esc)$ks[2],ks.boot(w.p,w.esc)$ks[2]))
	cat("PV for KS-test:\n",ks.ck.pv[1],"(P vs CSC)\n",ks.ck.pv[2],"(CSC vs ESC)\n",ks.ck.pv[3],"(P vs ESC)\n\n")
	t.ck.pv = c(t.test(w.p,w.csc)$p.value,t.test(w.csc,w.esc)$p.value,t.test(w.p,w.esc)$p.value)
	cat("PV for T-test:\n",t.ck.pv[1],"(P vs CSC)\n",t.ck.pv[2],"(CSC vs ESC)\n",t.ck.pv[3],"(P vs ESC)\n\n")
	q.out = cbind(q.out,ks.ck.pv)
}
q.out = q.out[,2:ncol(q.out)]
colnames(q.out) = n.label
ck = q.out[1,]
ck.out = q.out[,ck<0.05]
colnames(ck.out)
rownames(q.out) = q.type

write.table(q.out,'test_g_core.txt',sep="\t", col.names =T, row.names=T, quote=F)

# ----- plot per gene ---- *
plot(w.p,col = w.col[i.p],pch = 1, ylim = c(min(w.ck)-1,max(w.ck)+5),xlim =c(1,34), main=sprintf("Gene expression - %s",n.label[gind]),xlab="",ylab="dynamic range")
lines(w.p,col='green')
points(w.csc,col = w.col[i.csc],pch = 2)
lines(w.csc,col='red')
points(w.esc,col = w.col[i.esc],pch = 3)
lines(w.esc,col='orange')

pvv = ks.ck.pv
txt.pCsc = sprintf("P vs T: pv = %.3f", pvv[1])
txt.cscEsc = sprintf("T vs E: pv = %.3f", pvv[2])
txt.pEsc = sprintf("P vs E: pv = %.3f", pvv[3])
txt = c(txt.pCsc,txt.cscEsc,txt.pEsc)
legend(1,max(w.ck)+5,txt,cex=0.65,text.col = 'black', bg ='white',bty="n")
legend(1.2,max(w.ck)+3.5,c("PTCs","TSLCs","ESCs"),cex=0.65,col=c('green','red','orange'),lty=1,text.col = 'black', bg ='white',bty="n")
legend(10,max(w.ck)+5,sample,cex = 0.38,col = w.col,pch = c(rep(1,length(i.p)),rep(2,length(i.csc)),rep(3,length(i.esc))),ncol = 5,text.col='black',bg ='white',bty="n")

