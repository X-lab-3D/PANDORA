args <- commandArgs(trailingOnly = TRUE)

path = args[1]
ligand = args[2]
allele = args[3]
pdfname = args[4]

aff_thr = 0

rdata = read.table(path,head=T)

data_top1 = subset(rdata, rdata$aff1>aff_thr)[,c(1,2,3,4)]

pos_1 = data_top1[,3]

minoff = min(c(0,pos_1))
maxoff = max(c(nchar(ligand)-7,pos_1))

mean_aff_1 <- tapply(data_top1$aff1,data_top1$start1,mean)

##calculate histogram
bk = seq(from=minoff-0.5,to=maxoff+0.5,by=1)
h1 <- hist(pos_1,breaks=bk,plot=FALSE)
most_votes = which.max(h1$counts)+minoff

##labels
labs=vector(length=maxoff-minoff+1)

i=minoff
lind=1
while (i<0) {
    labs[lind] = 'X'
    lind = lind+1
    i = i+1
}
for (i in 0:maxoff) {
   labs[lind] = substring(ligand,i+1,i+1)
   lind = lind+1
}

###highlight core in upper case
prev = tolower(substring(ligand, 1, most_votes-1))
core = toupper(substring(ligand, most_votes, most_votes+8))
aft =  tolower(substring(ligand, most_votes+9))
c_ligand = paste(prev,core,aft,sep="")

##plot histogram###
pdf(pdfname,width=8,height=8)

plot(h1, freq=FALSE, xlab="P1 position in core",xaxt="n", ylab = "Core reliability score", ylim=c(0, 1.1*max(h1$density)), labels=labs, main=paste(allele,c_ligand,sep = " - "), col="whitesmoke")
axis(1, at=h1$mids, seq(from=minoff,to=maxoff,by=1))

dev.off()