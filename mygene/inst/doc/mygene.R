### R code from vignette source 'mygene.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: mygene.Rnw:31-32
###################################################
library(mygene)


###################################################
### code chunk number 3: mygene.Rnw:35-36
###################################################
mg.getgene("1017", fields=c("name","symbol","taxid","entrezgene"))


###################################################
### code chunk number 4: mygene.Rnw:46-47
###################################################
mg.getgenes(c("1017","1018","ENSG00000148795"), species="human")


###################################################
### code chunk number 5: mygene.Rnw:61-62
###################################################
mg.query("cdk2", size=5)


###################################################
### code chunk number 6: mygene.Rnw:67-68
###################################################
mg.query("symbol:cdk2", species="human")


###################################################
### code chunk number 7: mygene.Rnw:79-80
###################################################
mg.querymany(c("1017", "695"), scopes="entrezgene", species="human")


###################################################
### code chunk number 8: mygene.Rnw:86-87
###################################################
mg.querymany(c("1017","695"), scopes="entrezgene", species="9606", returnall=TRUE)


###################################################
### code chunk number 9: mygene.Rnw:104-114
###################################################
xli<-c('DDX26B',
 'CCDC83',
 'MAST3',
 'FLOT1',
 'RPL11',
 'ZDHHC20',
 'LUC7L3',
 'SNORD49A',
 'CTSH',
 'ACOT8')


###################################################
### code chunk number 10: mygene.Rnw:119-120
###################################################
mg.querymany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 11: mygene.Rnw:127-130
###################################################
out<-mg.querymany(xli, scopes="symbol", fields="ensembl.gene", species="human")
out
out$ensembl.gene[[4]]


###################################################
### code chunk number 12: mygene.Rnw:138-145
###################################################
xli<-c('DDX26B',
 'CCDC83',
 'MAST3',
 'FLOT1',
 'RPL11',
 'Gm10494')
mg.querymany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 13: mygene.Rnw:150-159
###################################################
xli<-c('DDX26B',
 'CCDC83',
 'MAST3',
 'FLOT1',
 'RPL11',
 'Gm10494',
 '1007_s_at',
 'AK125780')



###################################################
### code chunk number 14: mygene.Rnw:164-168
###################################################
out<-mg.querymany(xli, scopes=c("symbol","reporter","accession"), 
             fields=c("entrezgene","uniprot"), species="human")
out
out$uniprot[[5]]


###################################################
### code chunk number 15: mygene.Rnw:179-181
###################################################
mg.querymany(xli, scopes=c("symbol","reporter","accession"), 
             fields=c("entrezgene","uniprot"), species="human", returnall=TRUE)


