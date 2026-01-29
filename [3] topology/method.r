#=========================================================================
# Yau-Hua Yu: function file
#
# METHOD FILE
# ------------------------------------------------------------------------
# 1. get.parent/child
# 2. get.nb
# 3. is.tree2,  find.tree2
# 4. find.connect, find.cp, find.cp.all (net.m,label,nbhood)
# 5. focal.pth.cp(ind,net.m,label,nbhood,ck.s,out), is.cut.n, cut.dist.cp
# 6. node.attr.foc
# 7. k.b.dis, k.b.wt.dis2, k.b.wt.dis.G, k.b.wt.dis.G.contra
# -------------------------------------------------------------------------

# -- load library,files
source("sub_func.r")

# -------------------------------------------------------------------------										# --- 1.--- #
# 1 
# -------------------------------------------------------------------------
get.parent <- function(in.node,m,label) {
	node = m[,in.node]
	node.out = label[node ==1]
	ind.self = grep(label[in.node],node.out)
	if (length(ind.self)>0) {
		node.out[ind.self]=NA
	}
	node.out=node.out[!is.na(node.out)]
	return(node.out)
}

get.parent.fix <- function(in.node,m,label) {
	node = m[,in.node]
	node.out = label[node ==1]
	ck.node = label[in.node]
	ind.s = rep(F,length(node.out))
	ind.s[node.out==ck.node] = T
	
	if (length(ind.s)>0) {
		node.out[ind.s]=NA
	}
	node.out=node.out[!is.na(node.out)]
	return(node.out)
}

get.children <- function(in.node,m,label) {
	node = m[in.node,]
	node.out = label[node ==1]
	ind.self = grep(label[in.node],node.out)
	if (length(ind.self)>0) {
		node.out[ind.self]=NA
	}
	node.out=node.out[!is.na(node.out)]
	return(node.out)
}

get.children.fix <- function(in.node,m,label) {
	node = m[in.node,]
	node.out = label[node ==1]
	ck.node = label[in.node]
	ind.s = rep(F,length(node.out))
	ind.s[node.out==ck.node] = T
	
	if (length(ind.s)>0) {
		node.out[ind.s]=NA
	}
	node.out=node.out[!is.na(node.out)]
	return(node.out)
}


# -------------------------------------------------------------------------										# --- 2.--- #
# 2 
# -------------------------------------------------------------------------
get.nb <- function(in.node,m,label) {
	node.O = m[in.node,]
	node.out = label[node.O ==1]
	node.I = m[,in.node]
	node.in = label[node.I ==1]
	nb = c(node.out,node.in)
	if (length(nb)>1) {
		nb = removeRepeat(nb)
	}
	self = sprintf("^%s\\b$", label[in.node])
	m = regexpr(self,nb)
	out.nb = nb[m<0]
	return(out.nb)
}

get.nb.dg1 <- function(ind,net.m,label,out.nbhood) {
	nb = get.nb(ind,net.m,label)
	nb.dg = out.nbhood[match(nb,label)]
	ans = nb[nb.dg==1]
	return(ans)
}

get.nb.dg2 <- function(ind,net.m,label,out.nbhood) {
	nb = get.nb(ind,net.m,label)
	nb.dg = out.nbhood[match(nb,label)]
	ans = nb[nb.dg==2]
	return(ans)
}


# -------------------------------------------------------------------------										# --- 3.--- #
# 3 
# -------------------------------------------------------------------------
is.tree2 <- function(ind,net.m,label,out.nbhood) {
	ans = F
	
	nb = get.nb(ind,net.m,label)
	nb.dg = out.nbhood[match(nb,label)]
	nb.lg = length(nb)
	
	if (length(nb)>0) {
		if(sum(nb.dg)<=4) {
			ans = T
		}
		if(min(nb.dg)<=2) {
			ans = T
		}
	}	
	return(ans)
}

find.tree2 <- function(ind,net.m,label,out.nbhood) {
	ck = is.tree2(ind,net.m,label,out.nbhood)
	out = c()
	if (ck) {
		nb2 = get.nb.dg2(ind,net.m,label,out.nbhood)
		peri.nb = get.nb.dg1(ind,net.m,label,out.nbhood)
		m.nb2 = match(nb2,label)
		ck2 = lapply(m.nb2,is.tree2,net.m,label,out.nbhood)
		nb2 = nb2[unlist(ck2)]     #remove nodes connecting hubs
		m.nb2 = m.nb2[unlist(ck2)]
		out = c(label[ind],nb2,peri.nb)
		
		while (sum(unlist(as.numeric(ck2)))>0) {
			n.nb2 = unlist(lapply(m.nb2,get.nb.dg2,net.m,label,out.nbhood))
			peri.nb2 = unlist(lapply(m.nb2,get.nb.dg1,net.m,label,out.nbhood))
			m.n.nb2 = match(n.nb2,label)
			mck.n.nb2 = match(n.nb2,out)
			n.nb2 = n.nb2[is.na(mck.n.nb2)]
			m.n.nb2 = m.n.nb2[is.na(mck.n.nb2)]
			
			if (length(n.nb2)==0) {
				return(out)
			}
			if (length(n.nb2)>0) {
				ori1 = match(label[ind],n.nb2,nomatch=0)
				ori2 = match(label[ind],peri.nb2,nomatch=0)
				ori = c(ori1,ori2)
				if (sum(ori)>0) {
					ck2 = 0
					out = removeRepeat(append(out,n.nb2,peri.nb2))
				}
				if (sum(ori)==0) {
					ck2 = lapply(m.n.nb2,is.tree2,net.m,label,out.nbhood)
					out = removeRepeat(c(out,n.nb2,peri.nb2))
				}
			}
		}
		return(out)
	}
} 


# -------------------------------------------------------------------------										# --- 4.--- #
# 4 
# -------------------------------------------------------------------------
find.connect <- function(ind,m,label,nbhood,ck.s,collect) {
	out = collect
	ck.s[ind]=0
	tree = find.tree2(ind,m,label,nbhood)
	m.tree = match(tree,label)
	out = append(out,tree)
	cat("\n\nSTART of search: \t",label[ind],"\tQueue length: ",length(ck.s[ck.s==0]),"\n")
		
	nb = get.nb(ind,m,label)
	nnb = removeRepeat(append(nb,tree))
	m.nb = match(nnb,label,nomatch=0)
	out = ck.mergeQ(nb,out,out)
	ck.nb = ck.s[m.nb]
				
	ct =1
	while (sum(ck.nb)>0) {
		n.n = nnb[ck.nb!=0]
		cat("IN the loop: ",length(n.n),"\t Current Q [10]: ",n.n[1:10],"\n")
	
		ck.nb = ck.nb[ck.nb!=0]
		ord = order(ck.nb,decreasing=T)
		
		nt = n.n[ord[1]]
		ck.nb[ord[1]]=0
		ind.nt = match(nt,label)
		nt.tree = find.tree2(ind.nt,m,label,nbhood)
		nt.nb = get.nb(ind.nt,m,label)
		nnt.nb = removeRepeat(append(nt.nb,nt.tree))
		m.add = match(nnt.nb,out)
		add = nnt.nb[is.na(m.add)]

		#cat("length - before: \t",length(out),"\tTree length: \t",length(nt.tree),"\n")
		out = ck.mergeQ(nnt.nb,out,out)
		#cat("length - after: \t",length(out),"\n")
		ck.s[ind.nt]=0
				
		x = match(add,label)
		ck.x = ck.s[x]
		ind.x = x[ck.x>0]
		ck.x = ck.x[ck.x>0]
							
		nnb = c(n.n,label[ind.x])
		ck.nb = append(ck.nb,ck.x)
		ct = ct+1
		cat("LOOP COUNT: ",ct,"\t","new nb in Queue: ",length(nb),"\n")
						
	}
	return(removeRepeat(out))
}

find.cp <-find.connect; rm(find.connect)

find.cp.all <- function(net.m,label,nbhood) {
	# ck list
	ck.s = nbhood		# search Q starting pts
	ck.s[nbhood==1]=0	# remove peri genes from search Q
	ck.z = ck.s
	ck.f = rep(T,length(label))
	o = order(ck.z,decreasing=T)
	out = list()
	i = 1
	while (sum(ck.z)>0) {
		o = order(ck.z,decreasing=T)
		nt = label[o[1]]
		ind.nt = match(nt,label)
		ck.z[o[1]]=0
		cp=find.cp(ind.nt,net.m,label,nbhood,ck.s,nt)
		out[[i]] = cp
		ind.cp = match(cp,label)
		ck.z[ind.cp]=0
		ck.f[ind.cp]=F
		i=i+1
		cat("find.cp.all Count:",i,"\tNode in the loop: ",nt,"  cp length: ",length(cp),"\n")
	}
	# ck orphan, isolates (n=2)
	lg = length(out)
	orphan = label[ck.f]
	orph.dg = nbhood[ck.f]
	orph.ind = c(1:length(label))[ck.f]
	ck.o = rep(1,length(orphan))
	s = 1
	while (sum(ck.o)>0) {
		orph.dg = orph.dg[ck.o==1]
		orphan = orphan[ck.o==1]
		orph.ind = orph.ind[ck.o==1]
		ck.o = ck.o[ck.o==1]
		or = order(orph.dg,decreasing=T)
		dg = orph.dg[or[1]]
		if (dg==0) {
			out[[lg+s]] = orphan[or[1]]
			ck.o[or[1]]=0
			s = s+1
		}
		if (dg==1) {
			orph.nb = get.nb(orph.ind[or[1]],net.m,label)
			out[[lg+s]] = append(orphan[or[1]],orph.nb)
			m.o.nb = match(orph.nb,orphan)
			ck.o[or[1]]=0
			ck.o[m.o.nb]=0
			s = s+1			
		}
	}
	lg = sapply(out,length)
	out$lg = lg
	return(out)
}


# -------------------------------------------------------------------------										# --- 5.a --- #
# 5.a 
# -------------------------------------------------------------------------
focal.pth.cp <- function(ind,net.m,label,nbhood,ck.s,out) {						
	ans = list()
	node = label[ind]
	#nb = get.nb(ind,net.m,label)
	cp.old = find.cp(ind,net.m,label,nbhood,ck.s,out)
	cp.size = length(cp.old)
	m.cp = match(cp.old,label)
	net.m = net.m[m.cp,m.cp]
	label = label[m.cp]
	ck.s = ck.s[m.cp]
	nbhood = nbhood[m.cp]
	
	ind = match(node,label)	
	lg = length(label)
	n.nb = get.nb(ind,net.m,label)
	xd = c()
	if (ind==1) {
		xd = c(2:lg)
	}
	if (ind == lg) {
		xd = c(1:(lg-1))
	}
	if ((ind!=1)&&(ind!= lg)) {
		xd = c(1:(ind-1),(ind+1):lg)
	}
	m2 = net.m[xd,xd]
	lab = label[xd]
	m2.cks = ck.s[xd]
	m2.nbhood = nbhood[xd]
	
	nb.walk.out = find.cp.all(m2,lab,m2.nbhood)
			
	# calculate stats
	ans$pth = 0
	ans$cp.n = 1
	m.lg = cp.size -1	# m.lg =>start number: total node in the component[cp.size] - this node[removed]
	out.lg = 0
	x = nb.walk.out
	x.size = nb.walk.out$lg
	ans$cp.size = x.size
	x.cp = x[1:(length(x)-1)]
	ans$cp = list(x.cp)
	lg.x.cp = length(x.cp)	# should be the same with length(x.size)
	if(lg.x.cp>1) {
		x.size = sort(x.size,decreasing=T)	#sorting
		ans$cp.n = lg.x.cp
		for(i in 1:lg.x.cp) {	        
			x.lg = x.size[i]  
			m.lg = m.lg - x.lg
			out.lg = x.lg*m.lg
			ans$pth = ans$pth+out.lg
		}
	}
	ans$path = ans$pth/lg.x.cp
	return(ans$path)
}

is.cut.n <- function(ind,net.m,label,nbhood,ck.s,out) {
	ans = T
	node = label[ind]
	nb = get.nb(ind,net.m,label)
	#--- ck 1: periphery -> return 0
	node.dg = nbhood[ind]
	if (node.dg==1) {
		ans = F
	}
	#--- ck 2: if it is a cut node ;; NOT --> return 0
	cp.old = find.cp(ind,net.m,label,nbhood,ck.s,out)
	m.cp = match(cp.old,label)
	net.m = net.m[m.cp,m.cp]
	label = label[m.cp]
	ck.s = ck.s[m.cp]
	nbhood = nbhood[m.cp]
	lg = length(label)	# OLD cp size
	if (min(nbhood)>=(lg-1)/2) {
		ans = F
	}
	return(ans)
}

cut.dist.cp <- function(ind,net.m,label,nbhood,ck.s,out) {						
	ans = 0
	node = label[ind]
	cp.old = find.cp(ind,net.m,label,nbhood,ck.s,out)
	m.cp = match(cp.old,label)
	net.m = net.m[m.cp,m.cp]
	label = label[m.cp]
	ck.s = ck.s[m.cp]
	nbhood = nbhood[m.cp]
	lg = length(label)	# OLD cp size
	ind = match(node,label)	# ind of cut-node in its block (component)
	n.nb = get.nb(ind,net.m,label)
		xd = c()
		if (ind==1) {
			xd = c(2:lg)
		}
		if (ind == lg) {
			xd = c(1:(lg-1))
		}
		if ((ind!=1)&&(ind!= lg)) {
			xd = c(1:(ind-1),(ind+1):lg)
		}
	m2 = net.m[xd,xd]
	lab = label[xd]
	m2.cks = ck.s[xd]
	m2.nbhood = nbhood[xd]
	if (sum(m2.cks)>0) {
		new.cp = find.cp.all(m2,lab,m2.nbhood)
				
		# calculate stats
		ans.pth = 0
			
		m.lg = lg -1	# m.lg =>start number: total node in the component[cp.size] - this node[removed]
		out.lg = 0
		new.cp.lg = length(new.cp$lg)
		new.cp.size = new.cp$lg
		# ans = new.cp
			
		x.size = sort(new.cp.size,decreasing=T)	#sorting
		for(i in 1:new.cp.lg) {	        
			x.lg = x.size[i]  
			m.lg = m.lg - x.lg
			out.lg = x.lg*m.lg
			ans.pth = ans.pth+out.lg
		}
		ans = ans.pth/new.cp.lg
	}
	return(ans)
}


# -------------------------------------------------------------------------										# --- 6.--- #
# 6 
# -------------------------------------------------------------------------
node.attr.foc <- function(rxn.a,rxn.b,label) {
	tm=rxn2mat(rxn.a,rxn.b,label)
	net.m =tm[[1]]
	
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
	# --- grouping of node
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
	focality.ind = match(foc.g,label)   		
	flow.ind = match(flow.g,label)					
	peri.ind = match(peri.g,label)						
	node.type[focality.ind] = 1
	node.t.name[focality.ind] = "focality"
	node.type[flow.ind] = 2
	node.t.name[flow.ind] = "flow"
	node.type[peri.ind] = 3
	node.t.name[peri.ind] = "peri"
	node.sum = cbind(label,node.type,node.t.name,out.nbhood,out.foc.cp)
	return(node.sum)
}


# -------------------------------------------------------------------------										# --- 7.--- #
# 7 
# -------------------------------------------------------------------------
k.b.dis <- function(probs) {
	probs = as.numeric(probs[drop=T])
	lg = length(probs)
	out = 0
	for(i in 1:lg) {
		a = probs[i]
		#a.w = wt[i]
		cat(i)
		for(j in 1:lg){
			b = probs[j]
			dis = a*log(a/b)
			#n.dis = dis*a.w
			out = out+dis
			cat("\t",j," ",dis,"\tout:",out,"\n")
		}
	}
	out = out/lg/lg		
	return(out)
}

k.b.wt.dis2 <- function(probs, wt) {
	probs = as.numeric(probs[drop=T])
	lg = length(probs)
	out = 0
	for(i in 1:lg) {
		a = probs[i]
		a.w = wt[i]
		cat(i)
		for(j in 1:lg){
			b = probs[j]
			dis = a*log(a/b)
			n.dis = dis*a.w
			out = out + n.dis
			cat("\t",j," ",dis,"\tout:",out,"\n")
		}
	}
	out = 2*out/lg/(lg-1)		
	return(out)
}

k.b.wt.dis.G <- function(probs,wt) 		{		# break down for each gene
	probs = as.numeric(probs[drop=T])
	lg = length(probs)
	out.ck = c()	
	
	for(i in 1:lg) {
		a = probs[i]
		a.w = wt[i]
		#cat(i," out=",out.ck,"\n") 
		out = 0
		for(j in 1:lg){
			b = probs[j]
			dis = a*log(a/b)
			n.dis = dis*a.w
			out = out + n.dis
			#cat("\t",j," ",dis,"\tout:",out,"\n")
		}
		out.ck = append(out.ck,out)
		cat("check => ",out.ck,"\n")
	}
	return(out.ck)    
}

k.b.wt.dis.G.contra <- function(probs,wt) {	 	# for symmetric kullback-leibler dist
	probs = as.numeric(probs[drop=T])
	lg = length(probs)
	out.ck = c()	
	
	for(i in 1:lg) {
		a = probs[i]
		out = 0
		for(j in 1:lg){
			b.w = wt[j]
			b = probs[j]
			dis = b*log(b/a)
			n.dis = dis*b.w
			out = out + n.dis
		}
		out.ck = append(out.ck,out)
		cat("check => ",out.ck,"\n")
	}
	return(out.ck)    
}



