### R code from vignette source 'mygene.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: mygene.Rnw:31-33
###################################################
source('~/Documents/Su_Lab/Su_Lab-clone/mygene/R/mygene.R')
source('~/Documents/Su_Lab/Su_Lab-clone/mygene/R/utils.R')


###################################################
### code chunk number 3: mygene.Rnw:36-37
###################################################
getGene("1017", fields=c("name","symbol","taxid","entrezgene"))


###################################################
### code chunk number 4: mygene.Rnw:47-48
###################################################
getGenes(c("1017","1018","ENSG00000148795"), species="human")


###################################################
### code chunk number 5: mygene.Rnw:62-63
###################################################
query(q="cdk2", size=5)


###################################################
### code chunk number 6: mygene.Rnw:73-74
###################################################
queryMany(c(1017, 695), scopes="entrezgene", species="human")


###################################################
### code chunk number 7: mygene.Rnw:89-99
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
### code chunk number 8: mygene.Rnw:104-105
###################################################
queryMany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 9: mygene.Rnw:112-115
###################################################
out<-queryMany(xli, scopes="symbol", fields="ensembl.gene", species="human")
out
out$ensembl.gene[[4]]


###################################################
### code chunk number 10: mygene.Rnw:123-130
###################################################
xli<-c('DDX26B',
 'CCDC83',
 'MAST3',
 'FLOT1',
 'RPL11',
 'Gm10494')
queryMany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 11: mygene.Rnw:135-144
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
### code chunk number 12: mygene.Rnw:149-153
###################################################
out<-queryMany(xli, scopes=c("symbol","reporter","accession"), 
             fields=c("entrezgene","uniprot"), species="human")
out
out$`uniprot.Swiss-Prot`[[5]]


###################################################
### code chunk number 13: mygene.Rnw:164-166
###################################################
queryMany(xli, scopes=c("symbol", "reporter", "accession"), 
          fields=c("entrezgene", "uniprot"), species='human', returnall=TRUE)


