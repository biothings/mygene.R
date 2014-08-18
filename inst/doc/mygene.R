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
### code chunk number 3: mygene.Rnw:35-41
###################################################
gene=getGene("1017", fields="all")
length(gene)
gene$name
gene$taxid
gene$uniprot
gene$refseq


###################################################
### code chunk number 4: mygene.Rnw:51-52
###################################################
getGenes(c("1017","1018","ENSG00000148795"))


###################################################
### code chunk number 5: mygene.Rnw:66-67
###################################################
query(q="cdk2", size=5)


###################################################
### code chunk number 6: mygene.Rnw:70-71
###################################################
query(q="NM_013993")


###################################################
### code chunk number 7: mygene.Rnw:81-83
###################################################
queryMany(c('1053_at', '117_at', '121_at', '1255_g_at', '1294_at'),
          scopes="reporter", species="human")


###################################################
### code chunk number 8: mygene.Rnw:98-108
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
### code chunk number 9: mygene.Rnw:113-114
###################################################
queryMany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 10: mygene.Rnw:121-124
###################################################
out<-queryMany(xli, scopes="symbol", fields="ensembl.gene", species="human")
out
out$ensembl.gene[[4]]


###################################################
### code chunk number 11: mygene.Rnw:132-139
###################################################
xli<-c('DDX26B', 
       'CCDC83', 
       'MAST3', 
       'FLOT1', 
       'RPL11', 
       'Gm10494')
queryMany(xli, scopes="symbol", fields="entrezgene", species="human")


###################################################
### code chunk number 12: mygene.Rnw:144-153
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
### code chunk number 13: mygene.Rnw:158-162
###################################################
out<-queryMany(xli, scopes=c("symbol", "reporter","accession"), 
             fields=c("entrezgene","uniprot"), species="human")
out
out$`uniprot.Swiss-Prot`[[5]]


###################################################
### code chunk number 14: mygene.Rnw:173-175
###################################################
queryMany(xli, scopes=c("symbol", "reporter", "accession"), 
          fields=c("entrezgene", "uniprot"), species='human', returnall=TRUE)


