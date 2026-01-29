# compile 17 gene lists

source("sub_func.r")

# pre-process myc list (myc1+myc2)
# data = read.table("g_escSig_myc.txt",sep="\t",header=T,strip.white=T)
# id = as.character(data[,1][drop=T])
# sym = as.character(data[,2][drop=T])
# g = removeRepeat(sym)
# nu = data[match(g,sym),]
# write.table(nu,"g_escSig_myc.txt",sep="\t",quote=F,col.names=T)


# --------------------------------------------------------------------------------------------- #
# read in file list
# --------------------------------------------------------------------------------------------- #
f.list = read.delim("file_list.txt",sep="\t",header=F,strip.white=T)
f.list = as.character(f.list[,1][drop=T])
f.name = sprintf("%s.txt",f.list)

g.nR.d = c()
g.nR.s = c()
for (i in 1:length(f.list)) {
	data = read.table(f.name[i],sep="\t",header=T,strip.white=T)
	id = as.character(data[,1][drop=T])
	sym = as.character(data[,2][drop=T])
	cat("reading file: ",f.name[i],"\n")
	cat("check length:	",length(id),"\t",length(removeRepeat(id)),"\n")
	cat("current list length: ",length(g.nR.d),"\n")
	g.id = removeRepeat(id)
	g.sym = sym[match(g.id,id)]
	cat("----------> NA:\t",length(g.id[is.na(g.id)]),"\n")
	cat("----------> NA:\t",length(g.sym[is.na(g.sym)]),"\n")
	ck = match(g.id,g.nR.d,nomatch=0)
	g.nR.d = append(g.nR.d,g.id[ck==0])
	g.nR.s = append(g.nR.s,g.sym[ck==0])
	cat("total length: \t\t",length(g.nR.d),"\n\n")
}

lg = length(g.nR.d)
ind.m = matrix(0,nrow=lg,ncol=length(f.list))
for (i in 1:length(f.list)) {
	data = read.table(f.name[i],sep="\t",header=T,strip.white=T)
	id = as.character(data[,1][drop=T])
	sym = as.character(data[,2][drop=T])	
	ind.m[match(id,g.nR.d),i]=1
}

l.cnt = apply(ind.m,1,sum)
summary(l.cnt)					# min = 1, max = 9

colnames(ind.m) = f.list
out = cbind(g.nR.d,g.nR.s,l.cnt,ind.m)

write.table(out,"g_L_mat.txt",sep="\t",quote=F,col.names=T,row.names=F)




