# ===========================================================================================================================================
# precompiled 
# surv analysis looping functions
# 1. report.nw.surv
# 2. report.g.surv
# 3. report.g.cor.invT
#
# ------------------------------------------------------------------------------------------------------ -------- network all ---------------- # 
# 1																									 
# ----------------------------------------------------------------------------------------------------------------------------------------------
report.nw.surv <- function (epx,mag,spec,score,surv.mns,surv.st) {
	lg = length(epx)
	out.score = epx
	p.sort = summary(out.score)
	p1.gen = rep(0,lg)
	p1.gen[out.score<=p.sort[[2]]] = 1
	p1.gen[out.score>p.sort[[2]]] = 2
	p1.gen[out.score>p.sort[[3]]] = 3
	p1.gen[out.score>p.sort[[5]]] = 4
	
	out.score = mag
	p.sort = summary(out.score)														     
	p2.gen = rep(0,length(out.score))
	p2.gen[out.score<=p.sort[[2]]] = 1
	p2.gen[out.score>p.sort[[2]]] = 2
	p2.gen[out.score>p.sort[[3]]] = 3
	p2.gen[out.score>p.sort[[5]]] = 4
	
	out.score = spec
	p.sort = summary(out.score)														    
	p3.gen = rep(0,length(out.score))
	p3.gen[out.score<=p.sort[[2]]] = 1
	p3.gen[out.score>p.sort[[2]]] = 2
	p3.gen[out.score>p.sort[[3]]] = 3
	p3.gen[out.score>p.sort[[5]]] = 4
	
	out.score = score
	p.sort = summary(out.score)														     
	p4.gen = rep(0,length(out.score))
	p4.gen[out.score<=p.sort[[2]]] = 1
	p4.gen[out.score>p.sort[[2]]] = 2
	p4.gen[out.score>p.sort[[3]]] = 3
	p4.gen[out.score>p.sort[[5]]] = 4
																								# --------- network all --------------- #
	svt.data.1 = as.data.frame(as.list(cbind(surv.mns,surv.st,p1.gen)))
	svt.data.2 = as.data.frame(as.list(cbind(surv.mns,surv.st,p2.gen)))
	svt.data.3 = as.data.frame(as.list(cbind(surv.mns,surv.st,p3.gen)))
	svt.data.4 = as.data.frame(as.list(cbind(surv.mns,surv.st,p4.gen)))
		
	svt.g1 = survdiff(Surv(surv.mns,surv.st)~  p1.gen, data = svt.data.1)
	pv1 = 1-pchisq(svt.g1[[5]],3)
		
	svt.g2 = survdiff(Surv(surv.mns,surv.st)~  p2.gen, data = svt.data.2)
	pv2 = 1-pchisq(svt.g2[[5]],3)
		
	svt.g3 = survdiff(Surv(surv.mns,surv.st)~  p3.gen, data = svt.data.3)
	pv3 = 1-pchisq(svt.g3[[5]],3)
		
	svt.g4 = survdiff(Surv(surv.mns,surv.st)~  p4.gen, data = svt.data.4)
	pv4 = 1-pchisq(svt.g4[[5]],3)
		
	ans = c(pv1,pv2,pv3,pv4)
	return(ans)
}

# ------------------------------------------------------------------------------------------------------------------- per gene ---------------- # 
# 2																									 
# ----------------------------------------------------------------------------------------------------------------------------------------------
report.g.surv <- function (test.g,m.epx,m.mag,m.spec,m.score,surv.mns,surv.st,outF) {
	sgt.pv.1 = c()
	sgt.pv.2 = c()
	sgt.pv.3 = c()
	sgt.pv.4 = c()
	lg = length(test.g)
	
	for (t in 1:lg) {
		gind = match(test.g[t],rownames(m.epx))
		out.score = m.epx[gind,]
		p.sort = summary(out.score)
		p1.gen = rep(0,lg)
		p1.gen[out.score<=p.sort[[2]]] = 1
		p1.gen[out.score>p.sort[[2]]] = 2
		p1.gen[out.score>p.sort[[3]]] = 3
		p1.gen[out.score>p.sort[[5]]] = 4
		
		out.score = m.mag[gind,]
		p.sort = summary(out.score)														     
		p2.gen = rep(0,length(out.score))
		p2.gen[out.score<=p.sort[[2]]] = 1
		p2.gen[out.score>p.sort[[2]]] = 2
		p2.gen[out.score>p.sort[[3]]] = 3
		p2.gen[out.score>p.sort[[5]]] = 4
		
		out.score = m.spec[gind,]
		p.sort = summary(out.score)														    
		p3.gen = rep(0,length(out.score))
		p3.gen[out.score<=p.sort[[2]]] = 1
		p3.gen[out.score>p.sort[[2]]] = 2
		p3.gen[out.score>p.sort[[3]]] = 3
		p3.gen[out.score>p.sort[[5]]] = 4

		out.score = m.score[gind,]
		p.sort = summary(out.score)														     
		p4.gen = rep(0,length(out.score))
		p4.gen[out.score<=p.sort[[2]]] = 1
		p4.gen[out.score>p.sort[[2]]] = 2
		p4.gen[out.score>p.sort[[3]]] = 3
		p4.gen[out.score>p.sort[[5]]] = 4
																									# ------- per gene --------------- #		
		svt.data.1 = as.data.frame(as.list(cbind(surv.mns,surv.st,p1.gen)))
		svt.data.2 = as.data.frame(as.list(cbind(surv.mns,surv.st,p2.gen)))
		svt.data.3 = as.data.frame(as.list(cbind(surv.mns,surv.st,p3.gen)))
		svt.data.4 = as.data.frame(as.list(cbind(surv.mns,surv.st,p4.gen)))
			
		svt.g1 = survdiff(Surv(surv.mns,surv.st)~  p1.gen, data = svt.data.1)
		pv1 = 1-pchisq(svt.g1[[5]],3)
		sgt.pv.1 = append(sgt.pv.1,pv1)
		
		svt.g2 = survdiff(Surv(surv.mns,surv.st)~  p2.gen, data = svt.data.2)
		pv2 = 1-pchisq(svt.g2[[5]],3)
		sgt.pv.2 = append(sgt.pv.2,pv2)
		
		svt.g3 = survdiff(Surv(surv.mns,surv.st)~  p3.gen, data = svt.data.3)
		pv3 = 1-pchisq(svt.g3[[5]],3)
		sgt.pv.3 = append(sgt.pv.3,pv3)
		
		svt.g4 = survdiff(Surv(surv.mns,surv.st)~  p4.gen, data = svt.data.4)
		pv4 = 1-pchisq(svt.g4[[5]],3)
		sgt.pv.4 = append(sgt.pv.4,pv4)
		
		cat("survival prediction for gene: \t\t",test.g[t],"\n\n")
		
	}
	ck = rep(0,length(sgt.pv.1))
	ck1 = ck
	ck1[sgt.pv.1<0.05]=1
	ck2 = ck
	ck2[sgt.pv.2<0.05]=1
	ck3 = ck 
	ck3[sgt.pv.3<0.05]=1
	ck4 = ck 
	ck4[sgt.pv.4<0.05]=1
	ck.sum = ck1 + ck2 + ck3 + ck4

	sgt = cbind(test.g,sgt.pv.1,sgt.pv.2,sgt.pv.3,sgt.pv.4,ck1,ck2,ck3,ck4,ck.sum)
	colnames(sgt) = c('gene','pv-epx','pv-mag','pv-spec','pv-score','ck1','ck2','ck3','ck4','ck-sum')
	write.table(sgt,outF,sep="\t", col.names =T, row.names=F, quote=F)
	return(sgt)
}


# ------------------------------------------------------------------------------------------------------ nw - corr - 1/T (surv) --------------- #
# 3																									 
# -----------------------------------------------------------------------------------------------------------------------------------------------
report.g.cor.invT <- function (test.g,m.epx,m.mag,m.spec,m.score,surv.mns,outF) {

	ind.ex = rep(T,length(surv.mns))
	ind.ex[surv.mns==0] = F
	surv.mns = surv.mns[ind.ex]
	
	ot.pv.1 = c()
	ot.c.1 = c()
	ot.pv.2 = c()
	ot.c.2 = c()
	ot.pv.3 = c()
	ot.c.3 = c()
	ot.pv.4 = c()
	ot.c.4 = c()
	lg = length(test.g)
	
	for (t in 1:lg) {
		gind = match(test.g[t],rownames(m.epx))
		
		g1 = m.epx[gind,ind.ex]
		w1 = cor.test(g1,(surv.mns/12)^-1)
		ot.pv.1 = append(ot.pv.1,w1$p.value)
		ot.c.1 = append(ot.c.1,w1$estimate)
		
		g2 = m.mag[gind,ind.ex]
		w2 = cor.test(g2,(surv.mns/12)^-1)
		ot.pv.2 = append(ot.pv.2,w2$p.value)
		ot.c.2 = append(ot.c.2,w2$estimate)
		
		g3 = m.spec[gind,ind.ex]
		w3 = cor.test(g3,(surv.mns/12)^-1)
		ot.pv.3 = append(ot.pv.3,w3$p.value)
		ot.c.3 = append(ot.c.3,w3$estimate)
																									#  nw - corr - 1/T (surv) --------------- #
		g4 = m.score[gind,ind.ex]
		w4 = cor.test(g4,(surv.mns/12)^-1)
		ot.pv.4 = append(ot.pv.4,w4$p.value)
		ot.c.4 = append(ot.c.4,w4$estimate)
	}
	ck = rep(0,length(ot.pv.1))
	ck1 = ck
	ck1[ot.pv.1<0.05]=1
	ck2 = ck
	ck2[ot.pv.2<0.05]=1
	ck3 = ck 
	ck3[ot.pv.3<0.05]=1
	ck4 = ck 
	ck4[ot.pv.4<0.05]=1
	ck.sum = ck1 + ck2 + ck3 + ck4

	ot = cbind(test.g,ot.pv.1,ot.c.1,ot.pv.2,ot.c.2,ot.pv.3,ot.c.3,ot.pv.4,ot.c.4,ck1,ck2,ck3,ck4,ck.sum)
	colnames(ot) = c('gene','pv-epx','cor-epx','pv-mag','cor-mag','pv-spec','cor-spec','pv-score','cor-score','ck1','ck2','ck3','ck4','ck-sum')
	write.table(ot,outF,sep="\t", col.names =T, row.names=F, quote=F)
	return(ot)
}



