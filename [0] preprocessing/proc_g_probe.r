# ----------------------------------------------------------------------------------------------------------------------------------------- #
# * re-organize g.matrix g.list g.probe
# ----------------------------------------------------------------------------------------------------------------------------------------- #

source("sub_func.r")
library(limma)

# ----- annotation ----- # 										
ann.dat <-  read.delim("hgU133p2_short.txt", sep="\t",header=T)
cat(c("The annotation data contains ", nrow(ann.dat), " rows, and ", ncol(ann.dat), "columns\n"))

# ----- extract 
ann.probe = as.character(ann.dat[,1][drop = T])
g.id = as.character(ann.dat[,2][drop = T])
g.sym = as.character(ann.dat[,3][drop = T])
g.chr = as.character(ann.dat[,4][drop = T])

g.omim = as.character(ann.dat[,5][drop = T])
g.ann.tsp = as.character(ann.dat[,6][drop = T])
g.ensb = as.character(ann.dat[,7][drop = T])
g.name = as.character(ann.dat[,8][drop = T])

# ----- txt processing ----- #
txt.gid = strsplit(g.id,"///")
g.id.num = unlist(lapply(txt.gid,length))
length(g.id.num[g.id.num==2])

tot.p2id = sum(g.id.num)						# 58117
tot.pb = c()
for (i in 1:length(g.id.num)) {
	lg = g.id.num[i]
	tot.pb = append(tot.pb,rep(ann.probe[i],lg))
}

txt.gid = unlist(txt.gid)
txt.gid = trimWhiteSpace(txt.gid)
nr.gid = removeRepeat(txt.gid)
length(nr.gid)									# 21422
xx = as.numeric(nr.gid)
nr.gid = nr.gid[!is.na(xx)]
length(nr.gid)									# 21421 


# ------- retrieve probes --------- #
t.line = array(NA,length(nr.gid))

for(i in 1:length(nr.gid)) {
	pat = sprintf("^%s$", nr.gid[i])
	s = regexpr(pat,txt.gid)
	s.g = txt.gid[s>0]
	s.lg = length(s.g)
	s.p = tot.pb[s>0]
		
	cat(s.p,"\n",s.g,"\n\n")
		
	t.line[i] = paste(nr.gid[i],"\t",s.lg,"\t",paste(s.p,collapse=","))
}

# -- retrive official gene names -- #
tmp.gid = nr.gid
chr = read.delim("summary_gene.txt", sep="\t",header=F,strip.white=T)
chr.gid = as.character(chr[,1][drop=T])
chr.sym = as.character(chr[,2][drop=T])
chr.map = as.character(chr[,4][drop=T])
chr.name = as.character(chr[,7][drop=T])

t.ind = match(tmp.gid,chr.gid,nomatch=0)
nr.gid = nr.gid[t.ind!=0]						
tmp.chr = chr[t.ind,]						
#write.table(tmp.chr,'temp_chr.txt',sep = '\t',quote=F,col.names =T ,row.names=F)

# ** remove and output **
t.line = t.line[t.ind!=0]
#cat(t.line, file= "hg133p2_nrGid_pb.txt", sep="\n")


# ---- Annotation file processing ------- #
out.gid = gsub("///"," ",g.id)
out.ensb = gsub("///",",",g.ensb)
out.ann.tsp = gsub("(.*). // (.*)","\\1",g.ann.tsp)
out.ann.tsp = gsub("(.*). // (.*)","\\1",out.ann.tsp)
out.ann.tsp = gsub("(.*). // (.*)","\\1",out.ann.tsp)
n.ann.tsp = gsub("(.*) identifier using (.*) transcripts","\\2",out.ann.tsp)

out = cbind(ann.probe,g.id.num,out.gid,g.sym,g.chr,g.omim,n.ann.tsp,out.ensb,g.name)
write.table(out,'hg133p2_ann.txt',sep = '\t',quote=F,col.names =T ,row.names=F)

# ======================================================================================================= #
# * reading gene lists 
glst <-  read.delim("g_matrix.txt", sep="\t",header=T)
gid = trimWhiteSpace(as.character(glst[,1][drop=T]))
g.m = glst[,4:ncol(glst)]

# * read-in gid-2-pb file again
g2pb = read.delim("hg133p2_gid.txt", sep="\t",header=F,strip.white=T)

g2pb.gid = trimWhiteSpace(as.character(g2pb[,1][drop=T]))
g2pb.num = as.numeric(g2pb[,2][drop=T])
g2pb.sym = as.character(g2pb[,3][drop=T])
g2pb.chr = as.character(g2pb[,4][drop=T])
g2pb.name = as.character(g2pb[,5][drop=T])
g2pb.pbs = trimWhiteSpace(as.character(g2pb[,6][drop=T]))

ca.ind = match(gid,g2pb.gid,nomatch=0)
length(ca.ind[ca.ind==0])							
													
ca.g = g2pb[ca.ind,]
ca.tot = sum(ca.g[,2])								# 13588
ca.n = as.numeric(ca.g[,2])

# * output ca.g.probe file * #
ca.pbs = as.character(ca.g[,6][drop=T])
ca.pbs = unlist(strsplit(ca.pbs,","))
ca.pbs = trimWhiteSpace(ca.pbs)
length(ca.pbs)

x.pbs = removeRepeat(ca.pbs)
length(x.pbs)										# 13048
cat(x.pbs, file= "ca_probe.txt", sep="\n")

write.table(ca.g,'ca_ann.txt',sep = '\t',quote=F,col.names =T ,row.names=F)

