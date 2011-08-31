fc = read.table(file='FC_upstream_100_1.fa_1000_CAA_counts.xls',sep="\t",quote="",header=T)
tp = read.table(file='TP_jgi_upstream_100_1.fa_1000_CAA_counts.xls',sep="\t",quote="",header=T)
pt = read.table(file='PT_jgi_upstream_100_1.fa_1000_CAA_counts.xls',sep="\t",quote="",header=T)

motifs = as.vector(sort(unique(fc$motif)))

analize = function(tab,motifs) {
  res = c()
  for(i in 1:length(motifs)) {
    motif = motifs[i]
    real = tab[tab$motif == motif,]$real[1]
    rand.avg = mean(tab[tab$motif == motif,]$random)
    rand.sd = sd(tab[tab$motif == motif,]$random)
    o.e = real / rand.avg
    z = (real-rand.avg)/rand.sd
    res = rbind(res,c(motif,real,rand.avg,rand.sd,o.e,z))
  }
  rownames(res) = res[,1]
  res = res[,2:ncol(res)]
  colnames(res) = c('real','random','sd','o.e','z')
  data.frame(res)
}

fc = analize(fc,motifs)
pt = analize(pt,motifs)
tp = analize(tp,motifs)


matrix = cbind(CAA = c(as.numeric(as.vector(fc$z[1])),as.numeric(as.vector(tp$z[1])),as.numeric(as.vector(pt$z[1]))),
               CAACAA = c(as.numeric(as.vector(fc$z[2])),as.numeric(as.vector(tp$z[2])),as.numeric(as.vector(pt$z[2]))),
               CAACAACAA = c(as.numeric(as.vector(fc$z[3])),as.numeric(as.vector(tp$z[3])),as.numeric(as.vector(pt$z[3]))),
               CAACAACAACAA = c(as.numeric(as.vector(fc$z[4])),as.numeric(as.vector(tp$z[4])),as.numeric(as.vector(pt$z[4]))),
               CAACAACAACAACAA = c(as.numeric(as.vector(fc$z[5])),as.numeric(as.vector(tp$z[5])),as.numeric(as.vector(pt$z[5]))))

rownames(matrix) = c('FC','TP','PT')

FC = 27074
PT = 10024
TP = 11386

pdf(file='CAA_enrichment_Z.pdf',paper='a4r',width=11,height=14,pointsize=7)

par(bg="white")

barplot(matrix[c(1,2,3),c(5,4,3,2,1)],
        beside=T,horiz=T,
        legend=c('FC  ','TP  ','PT  '),names.arg=c('(CAA)5','(CAA)4','(CAA)3','(CAA)2','(CAA)1'), args.legend=list(cex = 1.5),
        main = "Motif Enrichment Z scores in Putative Promoters", font.main = 4, cex.main = 1.5,
        sub = "Z scores calculated against the distribution of 1000 shufflings", font.sub = 4,
        cex.names = 1.5)

par(bg="lightgrey")

matrix = cbind(CAA = c(as.numeric(as.vector(fc$real[1])) / FC * 100, as.numeric(as.vector(tp$real[1])) / TP * 100, as.numeric(as.vector(pt$real[1])) / PT * 100),
               CAACAA = c(as.numeric(as.vector(fc$real[2])) / FC * 100, as.numeric(as.vector(tp$real[2])) / TP * 100, as.numeric(as.vector(pt$real[2])) / PT * 100),
               CAACAACAA = c(as.numeric(as.vector(fc$real[3])) / FC * 100, as.numeric(as.vector(tp$real[3])) / TP * 100, as.numeric(as.vector(pt$real[3])) / PT * 100),
               CAACAACAACAA = c(as.numeric(as.vector(fc$real[4])) / FC * 100, as.numeric(as.vector(tp$real[4])) / TP * 100, as.numeric(as.vector(pt$real[4])) / PT * 100),
               CAACAACAACAACAA = c(as.numeric(as.vector(fc$real[5])) / FC * 100, as.numeric(as.vector(tp$real[5])) / TP * 100, as.numeric(as.vector(pt$real[5])) / PT * 100))

rownames(matrix) = c('FC','TP','PT')

matrix = round(matrix,digits=2)

matrix = cbind(matrix,rownames(matrix))
matrix = rbind(matrix,c('(CAA)1','(CAA)2','(CAA)3','(CAA)4','(CAA)5',' '))
matrix = matrix[c(4,1,2,3),c(6,1,2,3,4,5)]

legend(600,15,matrix,ncol=ncol(matrix),title='Percentage of Promoters')



matrix = cbind(CAA = c(as.numeric(as.vector(fc$real[1])) , as.numeric(as.vector(tp$real[1])) , as.numeric(as.vector(pt$real[1])) ),
               CAACAA = c(as.numeric(as.vector(fc$real[2])) , as.numeric(as.vector(tp$real[2])) , as.numeric(as.vector(pt$real[2])) ),
               CAACAACAA = c(as.numeric(as.vector(fc$real[3])) , as.numeric(as.vector(tp$real[3])) , as.numeric(as.vector(pt$real[3])) ),
               CAACAACAACAA = c(as.numeric(as.vector(fc$real[4])) , as.numeric(as.vector(tp$real[4])) , as.numeric(as.vector(pt$real[4])) ),
               CAACAACAACAACAA = c(as.numeric(as.vector(fc$real[5])) , as.numeric(as.vector(tp$real[5])) , as.numeric(as.vector(pt$real[5])) ))

rownames(matrix) = c('FC','TP','PT')

matrix = round(matrix)

matrix = cbind(matrix,rownames(matrix))
matrix = rbind(matrix,c('(CAA)1','(CAA)2','(CAA)3','(CAA)4','(CAA)5',' '))
matrix = matrix[c(4,1,2,3),c(6,1,2,3,4,5)]

legend(600,10,matrix,ncol=ncol(matrix),title='Number of Promoters')

dev.off()

prop.test(c(as.numeric(as.vector(fc[1,1])),as.numeric(as.vector(tp[1,1]))),c(FC,TP),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[2,1])),as.numeric(as.vector(tp[2,1]))),c(FC,TP),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[3,1])),as.numeric(as.vector(tp[3,1]))),c(FC,TP),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[4,1])),as.numeric(as.vector(tp[4,1]))),c(FC,TP),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[5,1])),as.numeric(as.vector(tp[5,1]))),c(FC,TP),alternative='g')$p.value

prop.test(c(as.numeric(as.vector(fc[1,1])),as.numeric(as.vector(pt[1,1]))),c(FC,PT),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[2,1])),as.numeric(as.vector(pt[2,1]))),c(FC,PT),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[3,1])),as.numeric(as.vector(pt[3,1]))),c(FC,PT),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[4,1])),as.numeric(as.vector(pt[4,1]))),c(FC,PT),alternative='g')$p.value
prop.test(c(as.numeric(as.vector(fc[5,1])),as.numeric(as.vector(pt[5,1]))),c(FC,PT),alternative='g')$p.value


