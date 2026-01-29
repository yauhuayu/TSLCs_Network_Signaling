# ========================================================================================================
# process tables (data) from IPA into a single table
# * keep:  ipa_id | in_gid | fold | location | family | drug | ipa_gid
# * reference to summary_gene to retrieve official g.sym, g.name, g.chr
# * final columns (10)
#   ipa_id | ipa_gid | g.sym | g.chr | g.name | in_gid | fold | location | family | drug 
# ========================================================================================================

# -- load library,files
source("sub_func.r")

files = sprintf("net%d.txt",1:2)

col = c("ipa_id","in_gid","sig_fold","cc_Loc","g_family","tx_drug","ipa_gid","network")
out = col

# -- Load IPA gene table


for (i in 1:length(files)) {
	
	cat("Reading IPA gene-net table ......\t",files[i],"\t")
	
	inf = files[i]
	input <- read.delim(inf, sep="\t",header=F)
	input = input[3:nrow(input),]				# remove first line of IPA intro, 2nd line; header
	cat("nrow of gene-net table ......\t",nrow(input),"\n")
	
	input = input[,c(1,4:9)]
	input = cbind(input,rep(i,nrow(input)))
	out = rbind(out,input)
}
out = out[2:nrow(out),]
colnames(out) = col

gid = as.character(out[,7])
gid[gid==""] = 0
gid = as.numeric(gid)
out[,7] = gid
outs = out[order(gid,decreasing=T),]

ipa.id = as.character(outs[,1][drop=T])
ipa.id = gsub(" ","",ipa.id)
ipa.id = gsub("\\(","_",ipa.id)
ipa.id = gsub("\\)","",ipa.id)
outs[,1] = ipa.id

nr.ipa.id = removeRepeat(ipa.id)
nw = as.character(outs[,8][drop=T])
nr.nw = rep(0,length(nr.ipa.id))

for(i in 1:length(nr.ipa.id)) {
	pat = sprintf("^%s\\b$", nr.ipa.id[i])
	m = regexpr(pat,ipa.id)
	m.nw = nw[m>0]
	nr.nw[i] = paste(m.nw,collapse=" ")
}
nr.ind = match(nr.ipa.id,ipa.id)
nr.out = outs[nr.ind,]
nr.out[,8] = nr.nw
nr.gid = as.numeric(nr.out[,7][drop=T])

# -----------------
# read in data
# -----------------
# * retrive official gene symbols
tmp.gid = nr.gid
chr = read.delim("summary_gene.txt", sep="\t",header=F,strip.white=T)
chr.gid = as.character(chr[,1][drop=T])
chr.sym = as.character(chr[,2][drop=T])
chr.map = as.character(chr[,4][drop=T])
chr.name = as.character(chr[,7][drop=T])

t.ind = match(tmp.gid,chr.gid,nomatch=0)
out.add= nr.out[t.ind!=0,]							
add.chr = chr[t.ind,c(2,4,7)]		
out.add = cbind(out.add,add.chr)

write.table(out.add,"lvcom_icord.txt",sep="\t",col.names =T ,row.names=F, quote=FALSE)

out.cpd = nr.out[t.ind==0,]

write.table(out.cpd,"lvcom_icord.txt",sep="\t",col.names =T ,row.names=F, quote=FALSE, append=T)





