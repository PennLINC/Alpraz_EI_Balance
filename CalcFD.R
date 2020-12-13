args<-commandArgs(TRUE)
if (length(args)==0){
	stop("Need input file!!")
}

df = read.table(args[1],header=TRUE,na.strings="n/a")
FD <- mean(df$framewise_displacement,na.rm=T)
write.table(FD,file= "",row.names=F,col.names=F)