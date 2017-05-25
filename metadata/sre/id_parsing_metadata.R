# making large metadata file for ID separation
# sre 25/5/17
#brachy 1,2,3,4, and 6 are single barcode axe files. They are parsed slightly different at first
#plate 5 has plate identifiers as 96 gel cut (BR05) and 96 XT cut (BRX5) per Niccy's notes
#

f1=read.table('brachy_01.axe',head=T)
names(f1)[1]='Barcode1'
f1$Plate='BR01'
f1$Well=paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep='')
f1$Axefile='brachy_01.axe'

f2=read.table('brachy_02.axe',head=T)
names(f2)[1]='Barcode1'
f2$Plate='BR02'
f2$Well=paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep='')
f2$Axefile='brachy_02.axe'

f3=read.table('brachy_03.axe',head=T)
names(f3)[1]='Barcode1'
f3$Plate='BR03'
f3$Well=paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep='')
f3$Axefile='brachy_03.axe'

f4=read.table('brachy_04.axe',head=T)
names(f4)[1]='Barcode1'
f4$Plate='BR04'
f4$Well=paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep='')
f4$Axefile='brachy_04.axe'

f5=read.table('brachy_05.axe',head=F)
names(f5)[1:3]=c('Barcode1','Barcode2','ID')
f5$Plate=c(rep('BR05',96),rep('BR05-XT',96))
f5$Well=rep(paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep=''),2)
f5$Axefile='brachy_05.axe'

f6=read.table('brachy_06.axe',head=T)
names(f6)[1]='Barcode1'
f6$Plate='B06'
f6$Well=paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep='')
f6$Axefile='brachy_06.axe'

f910=read.table('brachy_09-10.axe',head=F)
names(f910)[1:3]=c('Barcode1','Barcode2','ID')
f910$Plate=c(rep('BR09',96),rep('BR10',96))
f910$Well=rep(paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep=''),2)
f910$Axefile='brachy_9-10.axe'

f1114=read.table('brachy_11-14.axe',head=F)
names(f1114)[1:3]=c('Barcode1','Barcode2','ID')
f1114$Plate=c(rep('BR11',96),rep('BR12',96),rep('BR13',96),rep('BR14',96))
f1114$Well=rep(paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep=''),4)
f1114$Axefile='brachy_11-14.axe'

f1520=read.table('brachy_15-20.axe',head=F)
names(f1520)[1:3]=c('Barcode1','Barcode2','ID')
f1520$Plate=c(rep('BR15',96),rep('BR18',96),rep('BR19',96),rep('BR20',96))
f1520$Well=rep(paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep=''),4)
f1520$Axefile='brachy_15-20.axe'

f1617=read.table('brachy_16-17.axe',head=F)
names(f1617)[1:3]=c('Barcode1','Barcode2','ID')
f1617$Plate=c(rep('BR16',96),rep('BR17',96))
f1617$Well=rep(paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep=''),2)
f1617$Axefile='brachy_16-17.axe'

f2732=read.table('brachy_27-32.axe',head=F)
names(f2732)[1:3]=c('Barcode1','Barcode2','ID')
f2732$Plate=c(rep('BR27',96),rep('BR28',96),rep('BR29',96),rep('BR30',96),rep('BR31',96),rep('BR32',96))
f2732$Well=rep(paste(rep(c('A','B','C','D','E','F','G','H'),each=12),rep(c('01','02','03','04','05','06','07','08','09','10','11','12'),8),sep=''),6)
f2732$Axefile='brachy_27-32.axe'

library(gtools)
all=smartbind(f1,f2,f3,f4,f5,f6,f910,f1114,f1520,f1617,f2732)
all=all[,c('Barcode1','Barcode2','ID','Axefile','Plate','Well')]
all$SampleID=paste(all$Plate,all$Well,sep='')

#show that sampleIDs are uniqiue
length(all$SampleID)
length(unique(all$SampleID))

#write it out
write.table(all,'brachy_gbs_total_metadata.txt',sep='\t',row.names=F,quote=F)
