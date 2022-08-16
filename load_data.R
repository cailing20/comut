load('data/input.RData')
result.list<-readRDS('data/screen_results.rds')
sketch.gp<-htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th('co-mutation', title = 'co-mutation'),
      th('overall.comut.frac', title = 'co-mutation frequency across all samples regardless of cancer type'),
      th('concurrent.pv.sig', title = "number of cancer types with this gene pair concurrently mutated as assessed by Fisher's exact test, with nominal p-value < 0.05"),
      th('concurrent.padj.sig', title = "number of cancer types with this gene pair concurrently mutated as assessed by Fisher's exact test, with Benjamini-Hochberg procedure adjusted p-value < 0.05"),
      th('exclusive.pv.sig', title = "number of cancer types with this gene pair mutually exclusively mutated as assessed by Fisher's exact test, with nominal p-value < 0.05"),
      th('exclusive.padj.sig', title = "number of cancer types with this gene pair mutually exclusively  mutated as assessed by Fisher's exact test, with Benjamini-Hochberg procedure adjusted p-value < 0.05")
    )
  )
))
sketch.type<-htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th('CANCER_TYPE_DETAILED', title = 'detailed cancer types'),
      th('avg.total', title = 'averaged total cases belonging to the specific cancer type (the total cases vary by gene pair due to variable coverage of gene panels)'),
      th('concurrent.pv.sig', title = "number of gene pairs from this cancer type that were found to be concurrently mutated as assessed by Fisher's exact test, with nominal p-value < 0.05"),
      th('concurrent.padj.sig', title = "number of gene pairs from this cancer type that were found to be concurrently mutated as assessed by Fisher's exact test, with Benjamini-Hochberg procedure adjusted p-value < 0.05"),
      th('exclusive.pv.sig', title = "number of gene pairs from this cancer type that were found to be mutually exclusively mutated as assessed by Fisher's exact test, with nominal p-value < 0.05"),
      th('exclusive.padj.sig', title = "number of gene pairs from this cancer type that were found to be mutually exclusively  mutated as assessed by Fisher's exact test, with Benjamini-Hochberg procedure adjusted p-value < 0.05")
    )
  )
))
sketch.full<-htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th('co-mutation', title = 'gene pair assessed for co-mutation'),
      th('CANCER_TYPE_DETAILED', title = 'detailed cancer types'),
      th('comutant', title = 'number of samples with mutations in both genes'),
      th('g1_mutant', title = "number of samples with mutation in gene 1"),
      th('g2_mutant', title = "number of samples with mutation in gene 2"),
      th('total', title = "total number of samples from this cancer type that were included in panels that covered both genes"),
      th('co.frac',title="fraction of samples with co-mutations"),
      th('g1.frac',title="fraction of samples with mutations in gene 1"),
      th('g2.frac',title="fraction of samples with mutations in gene 2"),
      th('co.pv', title = "nominal p-value from Fisher's exact test for concurrent mutations"),
      th('co.adj', title = "p-value from Fisher's exact test for concurrent mutations adjusted for multiple comparisons by Benjamini-Hochberg procedures"),
      th('exclusive.pv', title = "nominal p-value from Fisher's exact test for mutually exclusive mutations"),
      th('exclusive.padj', title = "p-value from Fisher's exact test for mutually exclusive mutations adjusted for multiple comparisons by Benjamini-Hochberg procedures"),
      th('sig.level',title = "indicate whether the co-mutation is significant only by nominal p-value < 0.05 or remain significant after controlling for multiple comparison (p.adj < .05)"),
      th('occurrence',title='co-mutation pattern (concurrent or mutually exclusive)')
      
    )
  )
))
sketch.adhoc<-htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th('CANCER_TYPE_DETAILED', title = 'detailed cancer types'),
      th('comutant', title = 'number of samples with mutations in both genes'),
      th('g1_mutant', title = "number of samples with mutation in gene 1"),
      th('g2_mutant', title = "number of samples with mutation in gene 2"),
      th('total', title = "total number of samples from this cancer type that were included in panels that covered both genes"),
      th('co.frac',title="fraction of samples with co-mutations"),
      th('g1.frac',title="fraction of samples with mutations in gene 1"),
      th('g2.frac',title="fraction of samples with mutations in gene 2"),
      th('g1_mut_given_g2_mut',title="fraction of mutations in gene 1 given in samples with mutations in gene2"),
      th('g2_mut_given_g1_mut',title="fraction of mutations in gene 2 given in samples with mutations in gene1"),
      th('co.pv', title = "nominal p-value from Fisher's exact test for concurrent mutations"),
      th('exclusive.pv', title = "nominal p-value from Fisher's exact test for mutually exclusive mutations"),
      th('co.adj', title = "p-value from Fisher's exact test for concurrent mutations adjusted for multiple comparisons by Benjamini-Hochberg procedures"),
      th('exclusive.padj', title = "p-value from Fisher's exact test for mutually exclusive mutations adjusted for multiple comparisons by Benjamini-Hochberg procedures"),
      th('sig.level',title = "indicate whether the co-mutation is significant only by nominal p-value < 0.05 or remain significant after controlling for multiple comparison (p.adj < .05)"),
      th('occurrence',title='co-mutation pattern (concurrent or mutually exclusive)')
      
    )
  )
))
beautify.dt<-function(dt,pv.key=NULL){
  if(is.data.frame(dt) & !is.data.table(dt)) {df<-T;dt.rname<-rownames(dt);setDT(dt)} else df<-F
  num.col<-which(apply(as.matrix(dt),2,FUN=function(x) identical(is.na(x),is.na(as.numeric(x)))))
  int.col<-which(apply(as.matrix(dt),2,FUN=function(x) identical(as.numeric(x),as.numeric(as.integer(x)))))
  pv.cols<-num.col[grep(paste(c('pv','padj',pv.key),collapse ='|'),colnames(dt)[num.col])]
  for(pv.col in pv.cols) dt[[pv.col]]<-signif(dt[[pv.col]],2)
  for(other.col in setdiff(num.col,c(pv.cols,int.col))) dt[[other.col]]<-round(dt[[other.col]],4)
  if(df){dt<-as.data.frame(dt);rownames(dt)<-dt.rname}
  dt
}
test.comut<-function(g1,g2,include.CNA){
  g1.platforms<-unique(genie.combined$SEQ_ASSAY_ID[genie.combined$Hugo_Symbol==g1])
  g2.platforms<-unique(genie.combined$SEQ_ASSAY_ID[genie.combined$Hugo_Symbol==g2])
  common.platforms<-intersect(g2.platforms,g1.platforms)
  # add onto sample tables: patient, g1_mutation, g2_mutation, g1_CN, g2_CN, g1_fusion, g2_fusion
  shared.platform.ind<-which(m.data$SEQ_ASSAY_ID%in%common.platforms)
  if(grepl('A',g.ref[gene==g1]$mut.type) & include.CNA) g1_CN<-CNA.sub[match(m.data[shared.platform.ind]$SAMPLE_ID,rownames(CNA.sub)),g1]>0 else g1_CN<-F
  if(grepl('D',g.ref[gene==g1]$mut.type) & include.CNA) g1_CN<-CNA.sub[match(m.data[shared.platform.ind]$SAMPLE_ID,rownames(CNA.sub)),g1]<0 else g1_CN<-F
  if(grepl('A',g.ref[gene==g2]$mut.type) & include.CNA) g2_CN<-CNA.sub[match(m.data[shared.platform.ind]$SAMPLE_ID,rownames(CNA.sub)),g2]>0 else g2_CN<-F
  if(grepl('D',g.ref[gene==g2]$mut.type) & include.CNA) g2_CN<-CNA.sub[match(m.data[shared.platform.ind]$SAMPLE_ID,rownames(CNA.sub)),g2]<0 else g2_CN<-F
  
  if(!g.ref[gene==g1]$mut.type%in%c('A',NA)) g1_fusion<-m.data[shared.platform.ind]$SAMPLE_ID%in%fusion.sub$Tumor_Sample_Barcode[fusion.sub$Hugo_Symbol==g1] else g1_fusion<-F
  if(!g.ref[gene==g2]$mut.type%in%c('A',NA)) g2_fusion<-m.data[shared.platform.ind]$SAMPLE_ID%in%fusion.sub$Tumor_Sample_Barcode[fusion.sub$Hugo_Symbol==g2] else g2_fusion<-F
  mut.mat<-cbind(g1_mutation=m.data[shared.platform.ind]$SAMPLE_ID%in%mutations$Tumor_Sample_Barcode[which(mutations$Hugo_Symbol==g1)],
                 g2_mutation=m.data[shared.platform.ind]$SAMPLE_ID%in%mutations$Tumor_Sample_Barcode[which(mutations$Hugo_Symbol==g2)],
                 g1_CN,g2_CN,g1_fusion,g2_fusion)
  m.df<-data.table(m.data[shared.platform.ind,c('CANCER_TYPE','CANCER_TYPE_DETAILED'),with=F],g1.alternation=apply(mut.mat[,grep('g1',colnames(mut.mat))],1,any),g2.alteration=apply(mut.mat[,grep('g2',colnames(mut.mat))],1,any))
  by_cancer.df<-m.df[,list(comutant=length(which(g1.alternation&g2.alteration)),g1_mutant=length(which(g1.alternation)),g2_mutant=length(which(g2.alteration)),total=.N,
                           co.pv=tryCatch({fisher.test(table(g1.alternation,g2.alteration),alternative = 'greater')$p.value},error=function(e) -1.0),
                           exlusive.pv=tryCatch({fisher.test(table(g1.alternation,g2.alteration),alternative = 'less')$p.value},error=function(e) -1.0)),by=CANCER_TYPE_DETAILED]
  by_cancer.df[co.pv<0,co.pv:=NA][exlusive.pv<0,exlusive.pv:=NA]
  by_cancer.df[,co.frac:=comutant/total][,g1.frac:=g1_mutant/total][,g2.frac:=g2_mutant/total][,g1_mut_given_g2_mut:=comutant/g1_mutant][,g2_mut_given_g1_mut:=comutant/g2_mutant]
  by_cancer.df[,co.padj:=p.adjust(co.pv,'BH')][,exlusive.padj:=p.adjust(exlusive.pv,'BH')]
  by_cancer.df<-by_cancer.df[order(co.pv),c(grep("pv|padj",colnames(by_cancer.df),value = T,invert = T),grep("pv|padj",colnames(by_cancer.df),value = T)),with=F]
  by_cancer.df[,sig.level:=c('not significant','nominal pv < .05','p.adj < .05')[(1+(co.pv<.05|exlusive.pv<.05)+(co.padj<.05|exlusive.padj<.05))]]
  by_cancer.df[,occurrence:='']
  by_cancer.df[co.pv<.05,occurrence:='co-mutation'][exlusive.pv<.05,occurrence:='mutually-exclusive']
  beautify.dt(by_cancer.df)
}
plot.gene<-function(g1,g2,by_cancer.df){
  x <- seq(0,1,length.out=101)
  df <- expand.grid(x=x, y=x) #grid for colors
  df$expected<-df$x*df$y
  by_cancer.df[,`log10(total)`:=round(log10(total),4)]
  ggplot()+
    geom_tile(data =df, aes(x, y,fill = expected),alpha=.2)+
    geom_point(by_cancer.df[total>4],mapping=aes(x=g1.frac,y=g2.frac,size=`log10(total)`,label=CANCER_TYPE_DETAILED),alpha=0.3,pch=1)+
    geom_point(by_cancer.df[total>4],mapping=aes(x=g1.frac,y=g2.frac,color=co.frac,size=`log10(total)`,label=CANCER_TYPE_DETAILED),alpha=0.3,pch=19)+
    geom_point(by_cancer.df[which(total>4&sig.level!='not significant')],mapping=aes(x=g1.frac,y=g2.frac,shape=occurrence,alpha=sig.level))+
    scale_alpha_discrete(guide='none')+
    scale_shape_manual(values=c('co-mutation'=24,'mutually-exclusive'=25))+
    scale_size(range = c(2,10),breaks = c(1,2,3,4),labels = c('10^1','10^2','10^3','10^4'),guide = 'none')+
    scale_color_distiller(palette = 'Spectral',limits=c(0,1),guide='none')+scale_fill_distiller(palette = 'Spectral')+
    theme_classic()+theme(legend.position = 'bottom')+
    xlab(paste('Fraction with',g1,'alteration'))+ylab(paste('Fraction with',g2,'alteration'))+ labs(fill = "Fraction with alterations in both genes",size="Number of patients")+
    ggtitle(paste0("Analysis of ",g1," and ",g2, " alterations by ",length(by_cancer.df[total>4][['total']])," cancer types in ",comma(sum(by_cancer.df[total>4][['total']]))," patients from the AACR GENIE 11.0 data"))
  
}
