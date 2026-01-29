#=========================================================================
# Yau-Hua Yu: function file
#
# 1. removeRepeat(list)
# 2. ck2.same(a,b)       
#    ck2.diff(a,b)      
# 3. removeEmptyStr(list)
# 4. clap.rxn.fq(mat)
#    clap.rxn.data(mat)
# 5. ext.list.M(inL,mat)  
#    ext.list.M2(inL,mat) 	** exclusive 
# 6. get.sym(id,g.id,g.sym)
# 7. matrix2rxn(r.mat,label)
# 8. rxn2mat(a.list, b.list, label)
#=========================================================================


library(MASS)


#===============
# function -- 1																										* ---- 1. removeRepeat ---- *
#===============
removeRepeat<- function(list) {
	ls = length(list)
	index = rep(F, ls)
	
	for (i in (1:ls)) {
		x = match (list[i],list)
		if (!is.na(x)) {
			index[x] = T
		}
	}
	n.list = list[index]
	return (n.list)
}

outRepeat<- function(list) {
	ls = length(list)
	index = rep(T, ls)
	
	for (i in (1:ls)) {
		x = match (list[i],list)
		if (!is.na(x)) {
			index[x] = F
		}
	}
	n.list = list[index]
	return (n.list)
}



#===============
# function -- 2																									* ---- 2. ck2.same/diff ---- *
#===============
ck2.same <-function(a.list,b.list) {
	lg.a = length(a.list)
	lg.b = length(b.list)
	if (lg.a>lg.b) {
		m = match(b.list,a.list,nomatch=0)
		same = b.list[m!=0]
		return(same)
	}
	else if (lg.a<lg.b) {
		m = match(a.list,b.list,nomatch=0)
		same = a.list[m!=0]
		return(same)
	}
}


ck2.diff <-function(a.list,b.list) {
	lg.a = length(a.list)
	lg.b = length(b.list)
	m.a = match(a.list,b.list,nomatch=0)
	diff.a = a.list[m.a==0]
	m.b = match(b.list,a.list,nomatch=0)
	diff.b = b.list[m.b==0]
	out = c(diff.a,diff.b)
	return(out)
}


	
#===============
# function -- 3																									* ---- 3. removeEmptyStr ---- *
#===============
removeEmptyStr<- function(d) {
	ld = length(d)
	dd =rep(F, ld)
	str = sub('\\s+$', '', c(d[1:ld]), perl = TRUE)	

	for (i in 1:ld) {
 		m =match(str[i], c(""), nomatch = 0)
 		if (m == 0) {
 			dd[i] = T
 		}
	}

	d = d[dd] 
	return (d)
}

#===============
# function -- 4																									* ---- 4. clap.rxn.fq/data ---- *
#===============
clap.rxns.fq <- function(mat) {
	a = as.character(mat[,1][drop=T])
	b = as.character(mat[,2][drop=T])
	r = as.character(mat[,3][drop=T])
	
	rxn = sprintf("%s_%s",a,b)
	nr = removeRepeat(rxn)
	nr = sort(nr)
	out = c("n.out","n.end","freq")
	
	for (j in 1:length(nr)) {
		pat = sprintf("^%s\\b$", nr[j])
		m = regexpr(pat,rxn)
		m.rxn = rxn[m>0]
		lg.mr = length(m.rxn)
		n.a = unlist(strsplit(nr[j],"_"))[1]
		n.b = unlist(strsplit(nr[j],"_"))[2]
		add = c(n.a,n.b,lg.mr)
		out = rbind(out,add)
	}
	return(out)
}
		

clap.rxns.data <- function(mat) {
	a = as.character(mat[,1][drop=T])
	b = as.character(mat[,2][drop=T])
	if (ncol(mat)<3) {
		cat("Error: no data to collapse !\n")
	}
	if (ncol(mat)==3) {
		r = cbind(mat[,3],rep("",nrow(mat)))
		output = c("node.a","node.b","freq","data","NA")
	}
	if (ncol(mat)>3) {
		r = mat[,3:ncol(mat)]
		output = c("node.a","node.b","freq",colnames(mat)[3:ncol(mat)])
	}
	rxn = sprintf("%s_%s",a,b)
	nr = removeRepeat(rxn)
	nr = sort(nr)
	
	
	for (j in 1:length(nr)) {
		out = c()
		out1 = c()
		pat = sprintf("^%s\\b$", nr[j])
		m = regexpr(pat,rxn)
		m.rxn = rxn[m>0]
		lg.mr = length(m.rxn)
		m.r = r[m>0,]
		if (lg.mr==1) {
			n.a = unlist(strsplit(nr[j],"_"))[1]
			n.b = unlist(strsplit(nr[j],"_"))[2]
			add = c(n.a,n.b,lg.mr)
			out = rbind(out,add)
			out.r = m.r
			out1 = rbind(out1,out.r)
		}
		if (lg.mr>1) {
			n.a = unlist(strsplit(nr[j],"_"))[1]
			n.b = unlist(strsplit(nr[j],"_"))[2]
			add = c(n.a,n.b,lg.mr)
			out = rbind(out,add)
			out.r = m.r
			out.r.clap = apply(out.r,2,paste,collapse=" ")
			out1 = rbind(out1,out.r.clap)
		}
		op = cbind(out,out1)
		output = rbind(output,op)
	}
	return(output)
}




#===============
# function -- 5																									* ---- 5. ext.list.M/M2 ---- *
#===============
ext.list.M <- function(inL,mat) {
	a = as.character(mat[,1][drop=T])
	b = as.character(mat[,2][drop=T])
	r = as.character(mat[,3][drop=T])
	
	ind.o = rep(F,nrow(mat))
	for (j in 1:length(inL)) {
		pat = sprintf("^%s\\b$", inL[j])
		m = regexpr(pat,a)
		mm = regexpr(pat,b)
		ind.o[m>0] = T
		ind.o[mm>0] = T
	}
	out = mat[ind.o,]
	return(out)
}

ext.list.M2 <- function(inL,mat) {
	a = as.character(mat[,1][drop=T])
	b = as.character(mat[,2][drop=T])
	r = as.character(mat[,3][drop=T])
	
	ind.o = rep(F,nrow(mat))
	ind.1 = rep(0,nrow(mat))
	ind.2 = rep(0,nrow(mat))
	
	for (j in 1:length(inL)) {
		pat = sprintf("^%s\\b$", inL[j])
		m = regexpr(pat,a)
		mm = regexpr(pat,b)
		ind.1[m>0] = 1
		ind.2[mm>0] = 1
	}
	ind = ind.1+ind.2
	ind.o[ind==2] = T
	out = mat[ind.o,]
	return(out)
}
	

#===============
# function -- 6																									* ---- 6. get.sym ---- *
#===============
get.sym <- function(id,g.id,g.sym) {
	m = match(id,g.id)
	sym = g.sym[m]
	return(sym)
}	



#===============
# function -- 7																									* ---- 7. matrix2rxn ---- *
#===============
matrix2rxn <- function(r2.m,lab) {
	
	line = c()
	r2.m = data.matrix(r2.m)
	# store info
	self = c()
	self.num = 0
	arc = c()
	arc.num = 0
	edge = c()
	edge.num = 0
	
	for (p in 1:nrow(r2.m)) {
		# ---- ck: self ---- *
		if (r2.m[p,p]==1) {
			self.num = self.num +1
			self[[self.num]] = c(p,p)
		}
		# -- ck: NON-self -- *
		start = p
		if (p < nrow(r2.m)) {
			start = p+1
		}
		end = nrow(r2.m)
		for (q in start:end) {
			if ((r2.m[p,q]==1)&&(r2.m[q,p]==1 )){
				edge.num = edge.num +1
				edge[[edge.num]] = c(p,q)
			}
			else if ((r2.m[p,q]==1)&&(r2.m[q,p]==0)) {
				arc.num = arc.num +1
				arc[[arc.num]] = c(p,q)
			}
			else if ((r2.m[p,q]==0)&&(r2.m[q,p]==1)) {
				arc.num = arc.num +1
				arc[[arc.num]] = c(q,p)
			}
			else if ((r2.m[p,q]==0)&&(r2.m[q,p]==0)) {
				arc.num = arc.num
				edge.num = edge.num
			}
		}
	}
	
	# **** ouput **** #
	self.lg = length(self)
	if (self.lg >0) {
		self.line = array(NA,self.lg)
		for (y in 1:self.lg) {
			self.a = lab[self[[y]][1]]
			self.line[y] = paste(c("self-",self.a,self.a),collapse="\t")
		}
		line = self.line
	}
	
	arc.lg = length(arc)
	if (arc.lg >0) {
		arc.line = array(NA,arc.lg)
		for (y in 1:arc.lg) {
			arc.a = lab[arc[[y]][1]]
			arc.b = lab[arc[[y]][2]]
			arc.line[y] = paste(c("uni-",arc.a,arc.b),collapse="\t")
		}
		line = append(line, arc.line)
	}
		
	edge.lg = length(edge)
	if (edge.lg >0) {
		edge.line = array(NA,edge.lg)
		for (y in 1:edge.lg) {
			edge.item = lab[unlist(edge[[y]])]
			edge.line[y] = paste(c("bi-",edge.item),collapse="\t")
		}
		line = append(line,edge.line)
	}
	rxn = matrix(unlist(lapply(line,strsplit,"\t")),ncol=3,byrow=T)
	return(rxn)
}



# ==============													
# function -- 8																									* ---- 8. rxn2mat ---- *
# ==============
rxn2mat <- function(a.list, b.list, label) {
	out = list()
	
	a.id = as.character(a.list)
	b.id = as.character(b.list)
	
	node.num = length(label)
	
	# the matrix follow the order of the Label
	mat = matrix(0, ncol=node.num,nrow=node.num)
		
	for (i in 1:length(a.id)) {
		node.a = as.character(a.id[i])
		node.b = as.character(b.id[i])
		r.num = match(node.a,label)	# source node
		c.num = match(node.b,label)	# target node
		mat[r.num,c.num]=1
	}
	out$mat = mat
	out$node = label
	return(out)
}


