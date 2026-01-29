# ==========================================
# process txt file to get network rxn
# ==========================================


# ====------ load library,files -----==== #
source("sub_func.r")

col = c("Tissue","net.num","start","end")
out = col


# ==== read in netData ==== #
nw.chr <- read.delim("lvcom_icord.txt", sep="\t",header=T,strip.white=T)
nw.col = c('ipa_id','in_gid','sig_fold','cc_Loc','g_family','tx_drug','ipa_gid','network','g_sym','chr','g_name')
ipa.id = as.character(nw.chr[,1][drop=T])
ipa.gid = as.character(nw.chr[,7][drop=T])
ipa.nw = as.character(nw.chr[,8][drop=T])
g.sym = as.character(nw.chr[,9][drop=T])


# ====------ load IPA RXN file -----==== #

for (n in 1:2) {
	
	#n= 1

	cat("Reading IPA RXN-net table ......\t",n,"\t")
	
	inf = sprintf("lvcom_svg%d_out.txt",n)
	data <- read.delim(inf, sep="\t",header=F,strip.white=T)
	d.t = as.character(data[,1][drop=T])

	edge = data[d.t=="path:",2:ncol(data)]
	cat("number of edges : RXN-net table ......\t",nrow(edge),"\n")
	node = data[d.t!="path:",1:3]
	
	S.x = as.character(edge[,1][drop=T])
	S.x =  as.numeric(gsub("M", "", S.x))
	S.y = as.character(edge[,2][drop=T])
	S.y = as.numeric(S.y)
	T.x = as.character(edge[,3][drop=T])
	T.x =  as.numeric(gsub("L", "", T.x))
	T.y = as.character(edge[,4][drop=T])
	T.y = as.numeric(T.y)	
	
	n.name = as.character(node[,1][drop=T])
	n.name = n.name[1:(length(n.name)-2)]
	n.parse.x = as.character(node[,2][drop=T])
	n.parse.y = as.character(node[,3][drop=T])
	
	first =  gsub("x=", "", n.parse.x)
	sec =  gsub("y=", "", n.parse.y)
		
	N.x = as.numeric(first[1:(length(first)-2)])
	N.y = as.numeric(sec[1:(length(sec)-2)])
	
	#for (k in 1:length(N.x)) {				# correct the center for node name (too long)
	#	if (nchar(n.name[k])>10) {
	#		lg = nchar(n.name[k])
	#		move = 0.5*lg*10
	#		N.x[k] = N.x[k]+move
	#	}
	#}
	out.s = rep(NA,length(S.x))
	out.t = rep(NA,length(T.x))
	
	for (i in 1:length(S.x)) {
		
		ck.s = rep(NA,length(N.x))
		for (a in 1:length(N.x)) {
			ck.s[a] = ((N.x[a]-S.x[i])^2)+((N.y[a]-S.y[i])^2)
		}
		source = n.name[match(min(ck.s),ck.s)]
		
		ck.t = rep(NA,length(N.x))
		for (b in 1:length(N.x)) {
			ck.t[b] = (N.x[b]-T.x[i])^2+(N.y[b]-T.y[i])^2
		}
		target = n.name[match(min(ck.t),ck.t)]
		
		if (source!=target) {
			out.s[i] = source
			out.t[i] = target
		}
		if (source == target) {
			opt.1 = ck.s[order(ck.s)[2]]+min(ck.t)
			opt.2 = min(ck.s) + ck.t[order(ck.t)[2]]
			if (opt.1>opt.2) {
				out.s[i] = source
				out.t[i] = n.name[match(ck.t[order(ck.t)[2]],ck.t)]
			}
			if (opt.1<opt.2) {
				out.t[i] = target
				out.s[i] = n.name[match(ck.s[order(ck.s)[2]],ck.s)]
			}
		}
			
	}
	
	# another ck (the same edge cannot appear twice)  ;; i.e. some node is not marked
	ck = append(out.s,out.t)
	ck2 = ck
	ind.ck = match(n.name,ck)
	miss = n.name[is.na(ind.ck)]
		
	if (length(miss)>0) {
		mis.m = match(miss,n.name)
		for (y in mis.m) {
			ck.s2 = rep(NA,length(S.x))
			ck.t2 = rep(NA,length(T.x))
			for (r in 1:length(S.x)) {
				ck.s2[r] = ((N.x[y]-S.x[r])^2)+((N.y[y]-S.y[r])^2)
				ck.t2[r] = ((N.x[y]-T.x[r])^2)+((N.y[y]-T.y[r])^2)
			}
			cck = append(ck.s2,ck.t2)
			ck2[match(min(cck),cck)] = n.name[y]
		}
		ind.ck2 = match(n.name,ck2)
		miss = n.name[is.na(ind.ck2)]
	}
	out.s = ck2[1:length(S.x)]
	out.t = ck2[(length(S.x)+1):length(ck)]
	
	
	
	# another ck (the same edge cannot appear twice)  ;; deeper ck !!
	ck3 = rep(NA,length(S.x))
	for(u in 1:length(S.x)) {
		ck3[u] = sprintf("%s%s",out.s[u],out.t[u])
	}
	ck.y = outRepeat(ck3)
	if (length(ck.y)>0) {
		m.y = grep(ck.y,ck3)
		ck3.mat = array(NA,dim = c(length(m.y),3))
		for(t in 1:length(m.y)) {
			v = m.y[t]
			ck.s3 = rep(NA,length(N.x))
			ck.t3 = rep(NA,length(N.x))
			for (r in 1:length(N.x)) {
				ck.s3[r] = ((N.x[r]-S.x[v])^2)+((N.y[r]-S.y[v])^2)
				ck.t3[r] = ((N.x[r]-T.x[v])^2)+((N.y[r]-T.y[v])^2)
			}
			opt.1 = ck.s3[order(ck.s3)[2]]+min(ck.t3)
			opt.2 = min(ck.s3) + ck.t3[order(ck.t3)[2]]
			if (opt.1>opt.2) {
				ck3.mat[t,1] = opt.2
				ck3.mat[t,2] = n.name[match(ck.t3[order(ck.t3)[2]],ck.t3)]
				ck3.mat[t,3] = c("target")
			}
			if (opt.1<opt.2) {
				ck3.mat[t,1] = opt.1
				ck3.mat[t,2]  = n.name[match(ck.s3[order(ck.s3)[2]],ck.s3)]
				ck3.mat[t,3] = c("source")
			}
		}
		ck.out.dis = ck3.mat[,1]
		change = match(min(ck.out.dis),ck.out.dis)
		change.n = ck3.mat[change,3]
		if (change.n=="target") {
			out.t[m.y[change]] = ck3.mat[change,2]
		}
		if (change.n=="source") {
			out.s[m.y[change]] = ck3.mat[change,2]
		}
	}
	
	
	stage = "lvcom_ipa"
	net = n
		
	out.stage = rep(stage,length(out.s))
	out.net = rep(n,length(out.s))
	
	out.n =  cbind(out.stage,out.net,out.s,out.t)
	out = rbind(out,out.n)
	
	
}

out= out[2:nrow(out),]
colnames(out) = col

out.a = as.character(out[,3][drop=T])
out.a = gsub("\\(","_",out.a)
out.a = gsub("\\)","",out.a)
out.a = gsub("\\*","",out.a)
out.b = as.character(out[,4][drop=T])
out.b = gsub("\\(","_",out.b)
out.b = gsub("\\)","",out.b)
out.b = gsub("\\*","",out.b)

out.a.gid = ipa.gid[match(out.a, ipa.id)]
out.b.gid = ipa.gid[match(out.b, ipa.id)]
out.a.sym = g.sym[match(out.a, ipa.id)]
out.b.sym = g.sym[match(out.b, ipa.id)]

out = cbind(out,out.a.gid,out.b.gid,out.a.sym,out.b.sym)

# *************************** remove Repetitive links *************************** #
tst = sprintf("%s_%s",out.a.sym,out.b.sym)
db = as.character(out[,2][drop=T])
nr.rxn = removeRepeat(tst)
nw.out = c('ed-type','a.sym','b.sym','num-net','nw-net')

for (i in 1:length(nr.rxn)) {
	pat = nr.rxn[i]
	m = grep(pat,tst)
	p.db = paste(db[m],collapse="-")
	p.n.db = length(m)	
	
	p = unlist(strsplit(pat,"_"))
	p.a = p[1]
	p.b = p[2]
	p.ad = c("lvcom_ipa",p.a,p.b,p.n.db,p.db)
	nw.out = rbind(nw.out,p.ad)
}
	
nw.out = nw.out[2:nrow(nw.out),]
colnames(nw.out) = c('ed-type','a.sym','b.sym','num-net','nw-net')

a.gid = ipa.gid[match(nw.out[,2],g.sym)]
b.gid = ipa.gid[match(nw.out[,3],g.sym)]
nnw.out = cbind(nw.out[,1:3],a.gid,b.gid,nw.out[,4:5])

write.table(nnw.out,"lvcom_nw_icord.txt",sep="\t",col.names =T ,row.names=F, quote=FALSE)




# ---------------------------------------------------------------------------- # 
# process cpd-grs               NOT DONE for lvcom
# ---------------------------------------------------------------------------- #

## 
cpd = ipa.id[ipa.gid==0]
ind.ka = rep(1,length(out.a))
ind.kb = rep(1,length(out.b))

ad.out = colnames(out)
for (i in 1:length(cpd)) {
	pat = sprintf("^%s\\b$", cpd[i])
	m.a = regexpr(pat,out.a)
	ind.ka[m.a>0] = 0
	m.b = regexpr(pat,out.b)
	ind.kb[m.b>0] = 0
	
	nb.ab = out.b[m.a>0]
	nb.ab.gid = out.b.gid[m.a>0]
	add.a = nb.ab[nb.ab.gid!=0]
	
	nb.ba = out.a[m.b>0]
	nb.ba.gid = out.a.gid[m.b>0]
	add.b = nb.ba[nb.ba.gid!=0]
	
	add = removeRepeat(c(add.a,add.b))
	add.gid = ipa.gid[match(add,ipa.id)]
	add.nw = ipa.nw[match(add,ipa.id)]
	add.gsym = g.sym[match(add,ipa.id)]
	
	for (x in 1:length(add)) {
		adx = add[x]
		adx.nw = add.nw[x]
		adx.gid = add.gid[x]
		adx.sym = add.gsym[x]
		for (y in 1:length(add)) {
			if (x!=y) {
				ady = add[y]
				ady.nw = add.nw[y]
				ady.gid = add.gid[y]
				ady.sym = add.gsym[y]
				add.rxn = c("Lung",paste(c(adx.nw,cpd[i],ady.nw),collapse="-"),adx,ady,adx.gid,ady.gid,adx.sym,ady.sym)
				ad.out = rbind(ad.out,add.rxn)
			}
		}
	}
}

ad.out = ad.out[2:nrow(ad.out),]

# ck: 1) remove repeats
tst = sprintf("%s_%s",ad.out[,7],ad.out[,8])
nr.tst = removeRepeat(tst)
ad.nr = ad.out[match(nr.tst,tst),]

# ck: 2) remove one of the A-B and B-A 
tst = sprintf("%s_%s",ad.nr[,7],ad.nr[,8])
ad.a = as.character(ad.nr[,7][drop=T])
ad.b = as.character(ad.nr[,8][drop=T])
ad.lab = removeRepeat(c(ad.a,ad.b))
ad.mat = rxn2mat(ad.a,ad.b,ad.lab)[[1]]
ad.rxn = matrix2rxn(ad.mat,ad.lab)
ad.rxn = ad.rxn[2:length(ad.rxn)]
ad.rxn = gsub("\t","_",ad.rxn)
ad.f.out = ad.nr[match(ad.rxn,tst),]


# ---------------- 
ind.k = ind.ka+ind.kb
out.ncpd = out[ind.k==2,]

final = rbind(out.ncpd,ad.f.out)

write.table(final,"lung_nw_IPA.txt",sep="\t",col.names =T ,row.names=F, quote=FALSE)
