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
### code chunk number 3: mygene.Rnw:36-42
###################################################
gene=getGene("1017", fields="all")
length(gene)
gene$name
gene$taxid
gene$uniprot
gene$refseq


###################################################
### code chunk number 4: mygene.Rnw:52-53
###################################################
getGenes(c("1017","1018","ENSG00000148795"))


###################################################
### code chunk number 5: mygene.Rnw:67-68
###################################################
query(q="cdk2", size=5)


###################################################
### code chunk number 6: mygene.Rnw:71-72
###################################################
query(q="NM_013993")


###################################################
### code chunk number 7: mygene.Rnw:82-84
###################################################
queryMany(c('1053_at', '117_at', '121_at', '1255_g_at', '1294_at'),
          scopes="reporter", species="human")


###################################################
### code chunk number 8: mygene.Rnw:99-109
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
### code chunk number 9: mygene.Rnw:114-115
###################################################
queryMany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 10: mygene.Rnw:122-125
###################################################
out<-queryMany(xli, scopes="symbol", fields="ensembl.gene", species="human")
out
out$ensembl.gene[[4]]


###################################################
### code chunk number 11: mygene.Rnw:133-140
###################################################
xli<-c('DDX26B', 
       'CCDC83', 
       'MAST3', 
       'FLOT1', 
       'RPL11', 
       'Gm10494')
queryMany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 12: mygene.Rnw:145-154
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
### code chunk number 13: mygene.Rnw:159-163
###################################################
out<-queryMany(xli, scopes=c("symbol", "reporter","accession"), 
             fields=c("entrezgene","uniprot"), species="human")
out
out$`uniprot.Swiss-Prot`[[5]]


###################################################
### code chunk number 14: mygene.Rnw:174-176
###################################################
queryMany(xli, scopes=c("symbol", "reporter", "accession"), 
          fields=c("entrezgene", "uniprot"), species='human', returnall=TRUE)


