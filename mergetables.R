#Args: filex columnsfilex filey columnsfiley outputname


args = commandArgs(trailingOnly=TRUE)

filex <- read.delim(args[1],header=F)
filey <- read.delim(args[3],header=F)
colsx<-as.numeric(strsplit(args[2],split=",",fixed=TRUE)[[1]])
colsy<-as.numeric(strsplit(args[4],split=",",fixed=TRUE)[[1]]) 
filexy<-merge(filex,filey,by.x=colsx,by.y=colsy)

write.table(filexy,args[5],sep="\t",row.names=F, quote=F,col.names=F)
