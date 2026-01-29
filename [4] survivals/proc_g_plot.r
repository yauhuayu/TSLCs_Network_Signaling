# ----------------------------------------------------------------------------------------------------------------------------------------- #
# proc_g_plot.r
# ------------------------------------------------------------------------------------------------------------------------------------------
#
# plot survival
# 1. plot.surv
# 2. plot.surv.q.a
# 3. plot.surv.q.ab
# -------------------------------------------																	 # ** --- 1 --- ** #
plot.surv <- function (score,surv.mns,surv.st,xtitle) {
	p.sort = summary(score)														    									
	p.gen = rep(0,length(score))
	p.gen[score<=p.sort[[2]]] = 1
	p.gen[score>p.sort[[2]]] = 2
	p.gen[score>p.sort[[3]]] = 3
	p.gen[score>p.sort[[5]]] = 4
	svt.data = as.data.frame(cbind(surv.mns,surv.st,p.gen))
	svt.data = as.data.frame(as.list(svt.data))
	fit <- survfit(Surv(surv.mns,surv.st)~  p.gen, data = svt.data, type = 'kaplan-meier')
	w = survdiff(Surv(surv.mns,surv.st)~  p.gen, data = svt.data)
	pv = 1-pchisq(w[[5]],3)
	xpv = sprintf("p-value: %.3f (K-M)",pv)
	y = summary(coxph(Surv(surv.mns, surv.st) ~ p.gen))
	ypv = sprintf("p-value: %.3f (CoxPH)",y$coefficients[5])			# wald test
	# y$logtest[3] (log-rank test)

	plot(fit, col = 1:4,xlab='months',ylab='probability')
	title(main= xtitle,cex.main=1)
	legend(70,1.05, c("q1","q2","q3","q4",xpv,ypv),cex=0.5,col =1:4,lty = c(rep(1,4),0,0),text.col =c(1:4,1,1), bg = 'white',bty='n')
}

# -------------------------------------------																	 # ** --- 2 --- ** #
plot.surv.q.a <- function (score,surv.mns,surv.st,xtitle,a) {
	p.sort = summary(score)														    
	p.gen = rep(0,length(score))
	p.gen[score<=p.sort[[2]]] = 1
	p.gen[score>p.sort[[2]]] = 2
	p.gen[score>p.sort[[3]]] = 3
	p.gen[score>p.sort[[5]]] = 4
	
	p.gen2 = rep(0,length(score))
	p.gen2[p.gen==a] = 1
		
	svt.data = as.data.frame(cbind(surv.mns,surv.st,p.gen2))
	svt.data = as.data.frame(as.list(svt.data))
	fit <- survfit(Surv(surv.mns,surv.st)~  p.gen2, data = svt.data, type = 'kaplan-meier')
	w = survdiff(Surv(surv.mns,surv.st)~  p.gen2, data = svt.data)
	pv = 1-pchisq(w[[5]],1)
	xpv = sprintf("p-value: %.3f (K-M)",pv)
	
	plot(fit, col = c(2,4),lty = 2:3,xlab='years',ylab='survival')
	title(main=xtitle,cex.main=1)
	q = 1:4
	q.ind = rep(T,4)
	q.ind[q==a] = F
	q = q[q.ind]
	legend(5,0.4, c(sprintf("q(%d,%d,%d)",q[1],q[2],q[3]),sprintf("q%d",a),xpv), cex=0.75, col=c(2,4,1),lty =c(2,3,0),pch=c(3,3,-1),text.col = 'black', merge = TRUE, bg ='white',bty='n')
}

# -------------------------------------------																	 # ** --- 3 --- ** #
plot.surv.q.ab <- function (score,surv.mns,surv.st,xtitle,a,b) {
	p.sort = summary(score)														    
	p.gen = rep(0,length(score))
	p.gen[score<=p.sort[[2]]] = 1
	p.gen[score>p.sort[[2]]] = 2
	p.gen[score>p.sort[[3]]] = 3
	p.gen[score>p.sort[[5]]] = 4
	
	p.gen2 = rep(0,length(score))
	p.gen2[p.gen==a] = 1
	p.gen2[p.gen==b] = 1
		
	svt.data = as.data.frame(cbind(surv.mns,surv.st,p.gen2))
	svt.data = as.data.frame(as.list(svt.data))
	fit <- survfit(Surv(surv.mns,surv.st)~  p.gen2, data = svt.data, type = 'kaplan-meier')
	w = survdiff(Surv(surv.mns,surv.st)~  p.gen2, data = svt.data)
	pv = 1-pchisq(w[[5]],1)
		
	plot(fit, col = c(2,4),lty = 2:3,xlab='years',ylab='survival')
	title(main=xtitle,cex.main=1)
	q = 1:4
	q.ind = rep(T,4)
	q.ind[q==a] = F
	q.ind[q==b] = F
	q = q[q.ind]
	legend(5,0.4, c(sprintf("q%d,q%d",q[1],q[2]),sprintf("q%d,q%d",a,b),sprintf("p-value = %.3f",pv)), cex=0.75, col=c(2,4,1),lty =c(2,3,0),pch=c(3,3,-1),text.col = 'black', merge = TRUE, bg ='white')
}

