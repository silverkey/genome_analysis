t=read.table(file='CAACAACAA_KOG_2_FC_upstream_100_1_NSH_100_WL_9',sep="\t")
a = data.frame()
for(i in 1:nrow(t)) {
  r = t[i,]
  pval = prop.test(as.numeric(c(r[5],r[3])),as.numeric(c(r[4],r[2])),alternative='g')$p.value
  a = c(a,pval)
}
a = p.adjust(a,method='BH')

colnames(t) = c('class','total genes in KOG_2','total genes in class','(CAA)3 genes in KOG2','(CAA)3 genes in class')

t$adjusted.p.value = a

t[a<=0.05,]

selected = t[a<=0.05,]

selected = selected[order(selected$adjusted.p.value),]

write.table(selected,file='KOG_2_significant_CAA3.xls',sep="\t",quote=F,row.names=F)
