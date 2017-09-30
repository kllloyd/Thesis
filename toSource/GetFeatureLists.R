GetFeatureLists <- function(dataSetChoice,logHypStart,randomGeneList,chosenClinicalFeatures,covFuncForm,dimension){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Creates lists of features for use with Informed ARD models, with CuratedOvarian data
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	if(logHypStart[[1]]=='Inbuilt'&covFuncForm=='InformedARD'){
		switch(dataSetChoice,
				'OCGS'					={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,1)),'mean'=c(rep(0,dimension),0))},
				'OCGS+Clinical'			={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'SRGS'					={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,1)),'mean'=c(rep(0,dimension),0))},
				'SRGS+Clinical'			={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'OCGS+Rand'				={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'OCGS+Clin+Rand'		={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'SRGS+Rand'				={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'SRGS+Clin+Rand'		={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS'				={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS+Clin'		={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS+Rand'		={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS+Clin+Rand'	={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,4)),'mean'=c(rep(0,dimension),0))},
				'Random'				={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,1)),'mean'=c(rep(0,dimension),0))},
				'OCGSsublists'			={logHypStart <- list('noise'=log(0.2),'func'=log(0.8),'length'=log(rep(2,5)),'mean'=c(rep(0,dimension),0))})
	} else if(logHypStart[[1]]=='Inbuilt'&covFuncForm=='InformedARDV3'){
		switch(dataSetChoice,
				'OCGS'					={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,1)),'length'=log(rep(2,1)),'mean'=c(rep(0,dimension),0))},
				'OCGS+Clinical'			={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,2)),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'SRGS'					={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,1)),'length'=log(rep(2,1)),'mean'=c(rep(0,dimension),0))},
				'SRGS+Clinical'			={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,2)),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'OCGS+Rand'				={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,2)),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'OCGS+Clin+Rand'		={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,3)),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'SRGS+Rand'				={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,2)),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'SRGS+Clin+Rand'		={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,3)),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS'				={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,2)),'length'=log(rep(2,2)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS+Clin'		={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,3)),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS+Rand'		={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,3)),'length'=log(rep(2,3)),'mean'=c(rep(0,dimension),0))},
				'OCGS+SRGS+Clin+Rand'	={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,4)),'length'=log(rep(2,4)),'mean'=c(rep(0,dimension),0))},
				'Random'				={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,1)),'length'=log(rep(2,1)),'mean'=c(rep(0,dimension),0))},
				'OCGSsublists'			={logHypStart <- list('noise'=log(0.2),'func'=log(rep(0.8,5)),'length'=log(rep(2,5)),'mean'=c(rep(0,dimension),0))})
	}
	
	featureLists <- list('OCGS'		=c('AKT1','AKT2','AKT3','APAF1','BAD','BAX','BCL2','BCL2L1','BID','CFLAR','FAS','FASLG','HSPD1','HSPA1A','HSPA1L','HSP90AA1','HSP90AB1',
										'HSP90B1','BIRC2','IGF1','IGF1R','IGF2','IGF2R','IGFBP1','IGFBP2','DNAJC15','MCL1','MTOR','NFKB1','PIK3CA','PTEN','STAT3','BIRC5',
										'XIAP','ATP7B','ABCG2','CES1P1///CES1','CES2','NT5C2','DPYD','FPGS','H2AFX','GCLC','GCLM','GSTP1','SLC29A1','SLC29A2','ABCB1',
										'ABCC1','ABCC2','ABCC3','ABCC4','ABCC5','ABCC6','ABCC8','MVP','UMPS','RRM1','SOD1','TAP1','TAP2','ABCB4','TYMS','HPRT1','HMBS',
										'LINC00969///SDHA///SDHAP2///SDHAP1','TBP','ATM','BRCA1','ERCC1','ERCC2','MGMT','MLH1','MSH2','MSH6','RAD51','TOP1','TOP2A','TOP2B',
										'XPA','XRCC1','XRCC5','XRCC6','APC','TUBB2A///TUBB2B///TUBB3','PTGS2','EGFR','ERBB2','ERBB3','ERBB4','HIF1A','MKI67','CDKN2A',
										'CDKN1A','CDKN1B','TP53','VEGFA'),
						'SRGS'		=c('AGR2','MUTYH','AKAP12','TP53','TOP2A','FOXA2','SRC','SIVA1','ALDH9A1','LGR5','EHF','BAX','CES2','CPE','FGFBP1',
										'TUBB2A///TUBB2B///TUBB3///TUBB4A','ZNF12','RBM39','RFC3','GNPDA1','ANXA3','NFIB','ACTR3B','YWHAE','CYP51A1///LRRD1','HMGCS1',
										'ZMYND11','FADS2','SNX7','ARHGDIA','NDST1','DAP','ERCC8','GUCY1B3','HDAC1','HDAC2','IGFBP5','IL6','LSAMP','DGKZ','MYCBP','S100A10',
										'SLC1A3','NCOA1','TIAM1','VEGFA','RPL36','LBR','ABCB1','FASLG','TIMP1','FN1','TGFB1','XPA','POLH','ITGAE','ZNF200','COL3A1','ACKR3',
										'EPHB3','NBN','PCF11','DFNB31','BRCA2','AADAC','CD38','CHIT1','CXCR4','EFNB2','MECOM','FILIP1L','HSPB7','LRIG1','MMP1','PSAT1',
										'SDF2L1','TCF15','EPHB2','ETS1','TRIM27','MARK4','B4GALT5','ABCB10','AOC1'),
						'Clinical'	=c('grade','tumorstage','age'),
						'Random'	=randomGeneList,
						'OCGSapop' 	=c('AKT1','AKT2','AKT3','APAF1','BAD','BAX','BCL2','BCL2L1','BID','CFLAR','FAS','FASLG','HSPD1','HSPA1A','HSPA1L','HSP90AA1','HSP90AB1',
										'HSP90B1','BIRC2','IGF1','IGF1R','IGF2','IGF2R','IGFBP1','IGFBP2','DNAJC15','MCL1','MTOR','NFKB1','PIK3CA','PTEN','STAT3','BIRC5',
										'XIAP'),
						'OCGSpandd' =c('ATP7B','ABCG2','CES1P1///CES1','CES2','NT5C2','DPYD','FPGS','H2AFX','GCLC','GCLM','GSTP1','SLC29A1','SLC29A2','ABCB1',
										'ABCC1','ABCC2','ABCC3','ABCC4','ABCC5','ABCC6','ABCC8','MVP','UMPS','RRM1','SOD1','TAP1','TAP2','ABCB4','TYMS'),
						'OCGSdnare' =c('ATM','BRCA1','ERCC1','ERCC2','MGMT','MLH1','MSH2','MSH6','RAD51','TOP1','TOP2A','TOP2B','XPA','XRCC1','XRCC5','XRCC6'),
						'OCGSprlif' =c('APC','TUBB2A///TUBB2B///TUBB3','PTGS2','EGFR','ERBB2','ERBB3','ERBB4','HIF1A','MKI67','CDKN2A','CDKN1A','CDKN1B','TP53','VEGFA'),
						'OCGShkeep' =c('HPRT1','HMBS','LINC00969///SDHA///SDHAP2///SDHAP1','TBP'))

	switch(dataSetChoice,
		'OCGS'					={extraParam <- list('list1'=featureLists$OCGS)},
		'OCGS+Clinical'			={extraParam <- list('list1'=featureLists$OCGS,'list2'=featureLists$Clinical)},
		'SRGS'					={extraParam <- list('list1'=featureLists$SRGS)},
		'SRGS+Clinical'			={extraParam <- list('list1'=featureLists$SRGS,'list2'=featureLists$Clinical)},
		'OCGS+Rand'				={extraParam <- list('list1'=featureLists$OCGS,'list2'=featureLists$Random)},
		'OCGS+Clin+Rand'		={extraParam <- list('list1'=featureLists$OCGS,'list2'=featureLists$Clinical,'list3'=featureLists$Random)},
		'SRGS+Rand'				={extraParam <- list('list1'=featureLists$SRGS,'list2'=featureLists$Random)},
		'SRGS+Clin+Rand'		={extraParam <- list('list1'=featureLists$SRGS,'list2'=featureLists$Clinical,'list3'=featureLists$Random)},
		'OCGS+SRGS'				={extraParam <- list('list1'=featureLists$OCGS,'list2'=featureLists$SRGS)},
		'OCGS+SRGS+Clin'		={extraParam <- list('list1'=featureLists$OCGS,'list2'=featureLists$SRGS,'list3'=featureLists$Clinical)},
		'OCGS+SRGS+Rand'		={extraParam <- list('list1'=featureLists$OCGS,'list2'=featureLists$SRGS,'list3'=featureLists$Random)},
		'OCGS+SRGS+Clin+Rand'	={extraParam <- list('list1'=featureLists$OCGS,'list2'=featureLists$SRGS,'list3'=featureLists$Clinical,'list4'=featureLists$Random)},
		'Random'				={extraParam <- list('list1'=featureLists$Random)},
		'OCGSsublists' 			={extraParam <- list('list1'=featureLists$OCGSapop,'list2'=featureLists$OCGSpandd,'list3'=featureLists$OCGSdnare,'list4'=featureLists$OCGSprlif,'list5'=featureLists$OCGShkeep)})

	toReturn  <- list('logHypStart'=logHypStart,'extraParam'=extraParam,'dimension'=length(unlist(extraParam)))

	return(toReturn)

}