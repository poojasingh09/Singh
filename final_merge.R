path = "."
file_list <- dir(path, pattern ="5.h")

a <-read.table("DEXseq_GpCh_fullresults.sig0.05.h", header=TRUE, sep="\t")
b <-read.table("DEXseq_GpSd_fullresults.sig0.05.h", header=TRUE, sep="\t")
c <- merge(a,b,by="exon",all=TRUE)

d <-read.table("DEXseq_LdCh_fullresults.sig0.05.h", header=TRUE, sep="\t")
e <- merge(c,d,by="exon",all=TRUE)

f <-read.table("DEXseq_LdGp_fullresults.sig0.05.h", header=TRUE, sep="\t")
g <- merge(e,f,by="exon",all=TRUE)

h <-read.table("DEXseq_LdSd_fullresults.sig0.05.h", header=TRUE, sep="\t")
i <- merge(g,h,by="exon",all=TRUE)

j <-read.table("DEXseq_LdTd_fullresults.sig0.05.h", header=TRUE, sep="\t")
k <- merge(i,j,by="exon",all=TRUE)

l <-read.table("DEXseq_LdTm_fullresults.sig0.05.h", header=TRUE, sep="\t")
m <- merge(k,l,by="exon",all=TRUE)

n <-read.table("DEXseq_SdCh_fullresults.sig0.05.h", header=TRUE, sep="\t")
o <- merge(m,n,by="exon",all=TRUE)

p <-read.table("DEXseq_SdTm_fullresults.sig0.05.h", header=TRUE, sep="\t")
q <- merge(o,p,by="exon",all=TRUE)

r <-read.table("DEXseq_TdCh_fullresults.sig0.05.h", header=TRUE, sep="\t")
s <- merge(q,r,by="exon",all=TRUE)

t <-read.table("DEXseq_TdGp_fullresults.sig0.05.h", header=TRUE, sep="\t")
u <- merge(s,t,by="exon",all=TRUE)

v <-read.table("DEXseq_TdSd_fullresults.sig0.05.h", header=TRUE, sep="\t")
w <- merge(u,v,by="exon",all=TRUE)

x <-read.table("DEXseq_TdTm_fullresults.sig0.05.h", header=TRUE, sep="\t")
y <- merge(w,x,by="exon",all=TRUE)

z <-read.table("DEXseq_TmCh_fullresults.sig0.05.h", header=TRUE, sep="\t")
z1 <- merge(y,z,by="exon",all=TRUE)

z2 <-read.table("DEXseq_TmGp_fullresults.sig0.05.h", header=TRUE, sep="\t")
z3 <- merge(z1,z2,by="exon",all=TRUE)


write.table(z3, "merged.txt", sep="\t")
