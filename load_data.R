load('data/input.RData')
result.list<-readRDS('data/screen_results.rds')
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
plot.gene<-function(g1,g2,by_cancer.df,use.log=F){
  by_cancer.df[,`log10(total)`:=round(log10(total),3)]
  if(use.log) x <- 10^seq(log10(1e-3),log10(1),length.out=101) else x <- seq(0,1,length.out=101)
  df <- expand.grid(x=x, y=x) #grid for colors
  df$expected<-df$x*df$y
  by_cancer.df[,`log10(total)`:=round(log10(total),4)]
  g<-ggplot()+
    geom_tile(data =df, aes(x, y,fill = expected),alpha=.2)+
    geom_point(by_cancer.df[total>4],mapping=aes(x=g1.frac,y=g2.frac,size=`log10(total)`,label=CANCER_TYPE_DETAILED),alpha=0.3,pch=1)+
    geom_point(by_cancer.df[total>4],mapping=aes(x=g1.frac,y=g2.frac,color=co.frac,size=`log10(total)`,label=CANCER_TYPE_DETAILED),alpha=0.3,pch=19)+
    geom_point(by_cancer.df[which(total>4&sig.level!='not significant')],mapping=aes(x=g1.frac,y=g2.frac,shape=occurrence,alpha=sig.level))+
    scale_alpha_discrete(guide='none')+
    scale_shape_manual(values=c('co-mutation'=24,'mutually-exclusive'=25))+
    scale_size(range = c(2,10),breaks = c(1,2,3,4),labels = c('10^1','10^2','10^3','10^4'),guide = 'none')+
    theme_classic()+theme(legend.position = 'bottom')+
    xlab(paste('Fraction with',g1,'alteration'))+ylab(paste('Fraction with',g2,'alteration'))+ labs(fill = "Fraction with alterations in both genes",size="Number of patients")+
    ggtitle(paste0("Analysis of ",g1," and ",g2, " alterations by ",length(by_cancer.df[total>4][['total']])," cancer types in ",comma(sum(by_cancer.df[total>4][['total']]))," patients from the AACR GENIE 11.0 data"))
  if(use.log){
    g+scale_x_log10()+scale_y_log10()+scale_fill_distiller(palette = 'Spectral',breaks=c(0.001,0.01,0.1,1),trans='log')+scale_color_distiller(palette = 'Spectral',breaks=c(0.001,0.01,0.1,1),trans='log')
  }else{
    g+scale_color_distiller(palette = 'Spectral',limits=c(0,1),guide='none')+scale_fill_distiller(palette = 'Spectral')
  }
}
