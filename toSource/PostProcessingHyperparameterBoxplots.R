PostProcessingHyperparameterBoxplots <- function(unid,models,nReps,separateARD,plotHypRef){
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# K Lloyd 2016_08_09
	#--------------------------------------------------------------------------------------------------------------------------------------------#
	# Plot boxplots of final hyperparameter values for any number of repeats, for any number of GP or GPSurv models
	# If dimensionality if high, might want to plot ARD on new page, via separateARD
	# For synthetic data, can plot lines for generating hyperparameter values, via plotHypRef
	# Uses unids to load directly from files
	#--------------------------------------------------------------------------------------------------------------------------------------------#

	allModels <- c('GPNonSurvNoCens','GPSurvSqExp','GPSurvARD','GPSurvInformedARD','GPSurvInformedARDV2','GPSurvInformedARDV3','GPSurvInformedARDV4','GPSurvRF','GPSurvBIC','GPSurvInfMed','GPSurvInfUnif','GPSurvBICForward','GPSurvBICBackward')

	# Load required workspaces #
	if('GPNonSurvNoCens'%in%models) 	load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPNonSurvNoCens_Workspace.RData'))
	if('GPSurvSqExp'%in%models) 		load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvSqExp_Workspace.RData'))
	if('GPSurvARD'%in%models) 			load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvARD_Workspace.RData'))
	if('GPSurvInformedARD'%in%models) 	load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvInformedARD_Workspace.RData'))
	if('GPSurvInformedARDV2'%in%models) load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvInformedARDV2_Workspace.RData'))
	if('GPSurvInformedARDV3'%in%models) load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvInformedARDV3_Workspace.RData'))
	if('GPSurvInformedARDV4'%in%models) load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvInformedARDV4_Workspace.RData'))
	if('GPSurvRF'%in%models) 			load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvRF_Workspace.RData'))
	if('GPSurvBIC'%in%models) 			load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvBIC_Workspace.RData'))
	if('GPSurvInfMed'%in%models) 		load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvInfMed_Workspace.RData'))
	if('GPSurvInfUnif'%in%models) 		load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvInfUnif_Workspace.RData'))
	if('GPSurvBICForward'%in%models)	load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvBICForward_Workspace.RData'))
	if('GPSurvBICBackward'%in%models)	load(file=paste0(getwd(),'/Runs/',unid,'/',unid,'_outputStructureGPSurvBICBackward_Workspace.RData'))

	# Extract generating hyperparameters if relevant #
	if('GPNonSurvNoCens'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPNonSurvNoCens[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPNonSurvNoCens[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvSqExp'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvSqExp[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvSqExp[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvARD'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvARD[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvARD[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvInformedARD'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvInformedARD[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvInformedARD[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvInformedARDV2'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvInformedARDV2[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvInformedARDV2[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvInformedARDV3'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvInformedARDV3[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvInformedARDV3[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvInformedARDV4'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvInformedARDV4[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvInformedARDV4[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvRF'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvRF[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvRF[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvBIC'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvBIC[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvBIC[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvInfMed'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvInfMed[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvInfMed[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvInfUnif'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvInfUnif[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvInfUnif[[',1,']]$dataOptionsStructure$logHypGenerate')))
	} else if('GPSurvBICForward'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvBICForward[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvBICForward[[',1,']]$dataOptionsStructure$logHypGenerate')))
	}else if('GPSurvBICBackward'%in%models){
		dataSource <- eval(parse(text=paste0('outputStructureGPSurvBICBackward[[',1,']]$dataOptionsStructure$dataSource')))
		if(dataSource=='Generate') logHypGenerate <- eval(parse(text=paste0('outputStructureGPSurvBICBackward[[',1,']]$dataOptionsStructure$logHypGenerate')))
	}

	# Extract hyperparameters from each model #
	logHypChosen.GPNonSurvNoCens.noise 		<- numeric(nReps)
	logHypChosen.GPSurvSqExp.noise 			<- numeric(nReps)
	logHypChosen.GPSurvARD.noise			<- numeric(nReps)
	logHypChosen.GPSurvInformedARD.noise	<- numeric(nReps)
	logHypChosen.GPSurvInformedARDV2.noise	<- numeric(nReps)
	logHypChosen.GPSurvInformedARDV3.noise	<- numeric(nReps)
	logHypChosen.GPSurvInformedARDV4.noise	<- numeric(nReps)
	logHypChosen.GPSurvRF.noise 			<- numeric(nReps)
	logHypChosen.GPSurvBIC.noise 			<- numeric(nReps)
	logHypChosen.GPSurvInfMed.noise 		<- numeric(nReps)
	logHypChosen.GPSurvInfUnif.noise 		<- numeric(nReps)
	logHypChosen.GPSurvBICForward.noise 	<- numeric(nReps)
	logHypChosen.GPSurvBICBackward.noise 	<- numeric(nReps)

	if('GPNonSurvNoCens'%in%models)		logHypChosen.GPNonSurvNoCens.func 		<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPNonSurvNoCens[[',1,']]$logHypChosen$func')))))
	if('GPSurvSqExp'%in%models)			logHypChosen.GPSurvSqExp.func 			<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvSqExp[[',1,']]$logHypChosen$func')))))
	if('GPSurvARD'%in%models)			logHypChosen.GPSurvARD.func				<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvARD[[',1,']]$logHypChosen$func')))))
	if('GPSurvInformedARD'%in%models)	logHypChosen.GPSurvInformedARD.func		<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARD[[',1,']]$logHypChosen$func')))))
	if('GPSurvInformedARDV2'%in%models)	logHypChosen.GPSurvInformedARDV2.func	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARDV2[[',1,']]$logHypChosen$func')))))
	if('GPSurvInformedARDV3'%in%models)	logHypChosen.GPSurvInformedARDV3.func	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARDV3[[',1,']]$logHypChosen$func')))))
	if('GPSurvInformedARDV4'%in%models)	logHypChosen.GPSurvInformedARDV4.func	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARDV4[[',1,']]$logHypChosen$func')))))
	if('GPSurvRF'%in%models)			logHypChosen.GPSurvRF.func 				<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvRF[[',1,']]$logHypChosen$func')))))
	if('GPSurvBIC'%in%models)			logHypChosen.GPSurvBIC.func 			<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvBIC[[',1,']]$logHypChosen$func')))))
	if('GPSurvInfMed'%in%models)		logHypChosen.GPSurvInfMed.func 			<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInfMed[[',1,']]$logHypChosen$func')))))
	if('GPSurvInfUnif'%in%models)		logHypChosen.GPSurvInfUnif.func 		<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInfUnif[[',1,']]$logHypChosen$func')))))
	if('GPSurvBICForward'%in%models)	logHypChosen.GPSurvBICForward.func 		<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvBICForward[[',1,']]$logHypChosen$func')))))
	if('GPSurvBICBackward'%in%models)	logHypChosen.GPSurvBICBackward.func 	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvBICBackward[[',1,']]$logHypChosen$func')))))
	
	if('GPNonSurvNoCens'%in%models)		logHypChosen.GPNonSurvNoCens.length 	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPNonSurvNoCens[[',1,']]$logHypChosen$length')))))
	if('GPSurvSqExp'%in%models)			logHypChosen.GPSurvSqExp.length 		<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvSqExp[[',1,']]$logHypChosen$length')))))
	if('GPSurvARD'%in%models)			logHypChosen.GPSurvARD.length			<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvARD[[',1,']]$logHypChosen$length')))))
	if('GPSurvInformedARD'%in%models)	logHypChosen.GPSurvInformedARD.length	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARD[[',1,']]$logHypChosen$length')))))
	if('GPSurvInformedARDV2'%in%models)	logHypChosen.GPSurvInformedARDV2.length	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARDV2[[',1,']]$logHypChosen$length')))))
	if('GPSurvInformedARDV3'%in%models)	logHypChosen.GPSurvInformedARDV3.length	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARDV3[[',1,']]$logHypChosen$length')))))
	if('GPSurvInformedARDV4'%in%models)	logHypChosen.GPSurvInformedARDV4.length	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInformedARDV4[[',1,']]$logHypChosen$length')))))
	if('GPSurvRF'%in%models)			logHypChosen.GPSurvRF.length 			<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvRF[[',1,']]$logHypChosen$length')))))
	if('GPSurvBIC'%in%models)			logHypChosen.GPSurvBIC.length 			<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvBIC[[',1,']]$logHypChosen$length')))))
	if('GPSurvInfMed'%in%models)		logHypChosen.GPSurvInfMed.length 		<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInfMed[[',1,']]$logHypChosen$length')))))
	if('GPSurvInfUnif'%in%models)		logHypChosen.GPSurvInfUnif.length 		<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvInfUnif[[',1,']]$logHypChosen$length')))))
	if('GPSurvBICForward'%in%models)	logHypChosen.GPSurvBICForward.length 	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvBICForward[[',1,']]$logHypChosen$length')))))
	if('GPSurvBICBackward'%in%models)	logHypChosen.GPSurvBICBackward.length 	<- matrix(0,nrow=nReps,ncol=length(eval(parse(text=paste0('outputStructureGPSurvBICBackward[[',1,']]$logHypChosen$length')))))

	for(i in 1:nReps){
		if('GPNonSurvNoCens'%in%models){
			logHypChosen.GPNonSurvNoCens.noise[i]		<- eval(parse(text=paste0('outputStructureGPNonSurvNoCens[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPNonSurvNoCens.func[i]		<- eval(parse(text=paste0('outputStructureGPNonSurvNoCens[[',i,']]$logHypChosen$func')))
			logHypChosen.GPNonSurvNoCens.length[i,]		<- eval(parse(text=paste0('outputStructureGPNonSurvNoCens[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvSqExp'%in%models){
			logHypChosen.GPSurvSqExp.noise[i]			<- eval(parse(text=paste0('outputStructureGPSurvSqExp[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvSqExp.func[i]			<- eval(parse(text=paste0('outputStructureGPSurvSqExp[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvSqExp.length[i,]			<- eval(parse(text=paste0('outputStructureGPSurvSqExp[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvARD'%in%models){
			logHypChosen.GPSurvARD.noise[i]				<- eval(parse(text=paste0('outputStructureGPSurvARD[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvARD.func[i]				<- eval(parse(text=paste0('outputStructureGPSurvARD[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvARD.length[i,]			<- eval(parse(text=paste0('outputStructureGPSurvARD[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvInformedARD'%in%models){
			logHypChosen.GPSurvInformedARD.noise[i]		<- eval(parse(text=paste0('outputStructureGPSurvInformedARD[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvInformedARD.func[i]		<- eval(parse(text=paste0('outputStructureGPSurvInformedARD[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvInformedARD.length[i,]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARD[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvInformedARDV2'%in%models){
			logHypChosen.GPSurvInformedARDV2.noise[i]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV2[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvInformedARDV2.func[i]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV2[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvInformedARDV2.length[i,]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV2[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvInformedARDV3'%in%models){
			logHypChosen.GPSurvInformedARDV3.noise[i]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV3[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvInformedARDV3.func[i,]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV3[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvInformedARDV3.length[i,]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV3[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvInformedARDV4'%in%models){
			logHypChosen.GPSurvInformedARDV4.noise[i]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV4[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvInformedARDV4.func[i,]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV4[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvInformedARDV4.length[i,]	<- eval(parse(text=paste0('outputStructureGPSurvInformedARDV4[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvRF'%in%models){
			logHypChosen.GPSurvRF.noise[i]				<- eval(parse(text=paste0('outputStructureGPSurvRF[[',i,']]$logHypChosen$noise')))
		}
		if('GPSurvBIC'%in%models){
			logHypChosen.GPSurvBIC.noise[i]				<- eval(parse(text=paste0('outputStructureGPSurvBIC[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvBIC.func[i]				<- eval(parse(text=paste0('outputStructureGPSurvBIC[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvBIC.length[i,]			<- eval(parse(text=paste0('outputStructureGPSurvBIC[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvInfMed'%in%models){
			logHypChosen.GPSurvInfMed.noise[i]			<- eval(parse(text=paste0('outputStructureGPSurvInfMed[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvInfMed.func[i]			<- eval(parse(text=paste0('outputStructureGPSurvInfMed[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvInfMed.length[i,]		<- eval(parse(text=paste0('outputStructureGPSurvInfMed[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvInfUnif'%in%models){
			logHypChosen.GPSurvInfUnif.noise[i]			<- eval(parse(text=paste0('outputStructureGPSurvInfUnif[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvInfUnif.func[i]			<- eval(parse(text=paste0('outputStructureGPSurvInfUnif[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvInfUnif.length[i,]		<- eval(parse(text=paste0('outputStructureGPSurvInfUnif[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvBICForward'%in%models){
			logHypChosen.GPSurvBICForward.noise[i]		<- eval(parse(text=paste0('outputStructureGPSurvBICForward[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvBICForward.func[i]		<- eval(parse(text=paste0('outputStructureGPSurvBICForward[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvBICForward.length[i,]	<- eval(parse(text=paste0('outputStructureGPSurvBICForward[[',i,']]$logHypChosen$length')))
		}
		if('GPSurvBICBackward'%in%models){
			logHypChosen.GPSurvBICBackward.noise[i]		<- eval(parse(text=paste0('outputStructureGPSurvBICBackward[[',i,']]$logHypChosen$noise')))
			logHypChosen.GPSurvBICBackward.func[i]		<- eval(parse(text=paste0('outputStructureGPSurvBICBackward[[',i,']]$logHypChosen$func')))
			logHypChosen.GPSurvBICBackward.length[i,]	<- eval(parse(text=paste0('outputStructureGPSurvBICBackward[[',i,']]$logHypChosen$length')))
		}
	}

	# Construct vectors representing models to plot & labels #
	toBoxplot 		<- list()
	boxplotLabels 	<- list()
	if('GPNonSurvNoCens'%in%models){
		toBoxplot[[1]] 		<- c(1,1,1)
		boxplotLabels[[1]] 	<- c('noise','func','length')
	}
	if('GPSurvSqExp'%in%models){
		toBoxplot[[2]] 		<- c(2,2,2)
		boxplotLabels[[2]] 	<- c('noise','func','length')
	}
	if('GPSurvARD'%in%models&!separateARD){
		toBoxplot[[3]] 		<- c(3,3,rep(3,dim(logHypChosen.GPSurvARD.length)[2]))
		boxplotLabels[[3]] 	<- c('noise','func',paste0('length',1:dim(logHypChosen.GPSurvARD.length)[2]))
	}
	if('GPSurvInformedARD'%in%models){
		toBoxplot[[4]] 		<- c(4,4,rep(4,dim(logHypChosen.GPSurvInformedARD.length)[2]))
		boxplotLabels[[4]] 	<- c('noise','func',paste0('length',1:dim(logHypChosen.GPSurvInformedARD.length)[2]))
	}
	if('GPSurvInformedARDV2'%in%models){
		toBoxplot[[5]] 		<- c(5,5,rep(5,dim(logHypChosen.GPSurvInformedARDV2.length)[2]))
		boxplotLabels[[5]] 	<- c('noise','func',paste0('length',1:dim(logHypChosen.GPSurvInformedARDV2.length)[2]))
	}
	if('GPSurvInformedARDV3'%in%models){
		toBoxplot[[6]] 		<- c(6,rep(6,dim(logHypChosen.GPSurvInformedARDV3.func)[2]),rep(6,dim(logHypChosen.GPSurvInformedARDV3.length)[2]))
		boxplotLabels[[6]] 	<- c('noise',paste0('func',1:dim(logHypChosen.GPSurvInformedARDV3.length)[2]),paste0('length',1:dim(logHypChosen.GPSurvInformedARDV3.length)[2]))
	}
	if('GPSurvInformedARDV4'%in%models){
		toBoxplot[[7]] 		<- c(7,rep(7,dim(logHypChosen.GPSurvInformedARDV4.func)[2]),rep(7,dim(logHypChosen.GPSurvInformedARDV4.length)[2]))
		boxplotLabels[[7]] 	<- c('noise',paste0('func',1:dim(logHypChosen.GPSurvInformedARDV4.length)[2]),paste0('length',1:dim(logHypChosen.GPSurvInformedARDV4.length)[2]))
	}
	if('GPSurvRF'%in%models){
		toBoxplot[[8]] 		<- c(8)
		boxplotLabels[[8]] 	<- 'noise'
	}
	if('GPSurvBIC'%in%models){
		toBoxplot[[9]] 		<- c(9,9,9)
		boxplotLabels[[9]] 	<- c('noise','func','length')
	}
	if('GPSurvInfMed'%in%models){
		toBoxplot[[10]] 	<- c(10,10,10)
		boxplotLabels[[10]] <- c('noise','func','length')
	}
	if('GPSurvInfUnif'%in%models){
		toBoxplot[[11]] 	<- c(11,11,11)
		boxplotLabels[[11]] <- c('noise','func','length')
	}
	if('GPSurvBICForward'%in%models){
		toBoxplot[[12]] 		<- c(12,12,12)
		boxplotLabels[[12]] 	<- c('noise','func','length')
	}
	if('GPSurvBICBackward'%in%models){
		toBoxplot[[13]] 		<- c(13,13,13)
		boxplotLabels[[13]] 	<- c('noise','func','length')
	}
	toBoxplot 		<- unlist(toBoxplot)
	boxplotLabels 	<- unlist(boxplotLabels)

	# Compile hyperparameter data to plot #
	logHypChosenToBoxplot 	<- matrix(NA,nrow=nReps,ncol=length(toBoxplot))
	colNum 					<- 1
	if('GPNonSurvNoCens'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPNonSurvNoCens.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPNonSurvNoCens.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPNonSurvNoCens.length)[2])] 		<- logHypChosen.GPNonSurvNoCens.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPNonSurvNoCens.length)[2])
	}
	if('GPSurvSqExp'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvSqExp.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvSqExp.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvSqExp.length)[2])] 			<- logHypChosen.GPSurvSqExp.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvSqExp.length)[2])
	}
	if('GPSurvARD'%in%models&!separateARD){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvARD.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvARD.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvARD.length)[2])] 			<- logHypChosen.GPSurvARD.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvARD.length)[2])
	}
	if('GPSurvInformedARD'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvInformedARD.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvInformedARD.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvInformedARD.length)[2])] 	<- logHypChosen.GPSurvInformedARD.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvInformedARD.length)[2])
	}
	if('GPSurvInformedARDV2'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvInformedARDV2.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvInformedARDV2.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvInformedARDV2.length)[2])] 	<- logHypChosen.GPSurvInformedARDV2.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvInformedARDV2.length)[2])
	}
	if('GPSurvInformedARDV3'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvInformedARDV3.noise
		logHypChosenToBoxplot[,(colNum+1):(colNum+dim(logHypChosen.GPSurvInformedARDV3.func)[2])] 		<- logHypChosen.GPSurvInformedARDV3.func
		logHypChosenToBoxplot[,(colNum+1+dim(logHypChosen.GPSurvInformedARDV3.func)[2]):(colNum+dim(logHypChosen.GPSurvInformedARDV3.func)[2]+dim(logHypChosen.GPSurvInformedARDV3.length)[2])] 	<- logHypChosen.GPSurvInformedARDV3.length
		colNum 																							<- colNum+1+dim(logHypChosen.GPSurvInformedARDV3.func)[2]+dim(logHypChosen.GPSurvInformedARDV3.length)[2]
	}
	if('GPSurvInformedARDV4'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvInformedARDV4.noise
		logHypChosenToBoxplot[,(colNum+1):(colNum+dim(logHypChosen.GPSurvInformedARDV4.func)[2])] 		<- logHypChosen.GPSurvInformedARDV4.func
		logHypChosenToBoxplot[,(colNum+1+dim(logHypChosen.GPSurvInformedARDV4.func)[2]):(colNum+dim(logHypChosen.GPSurvInformedARDV4.func)[2]+dim(logHypChosen.GPSurvInformedARDV4.length)[2])] 	<- logHypChosen.GPSurvInformedARDV4.length
		colNum 																							<- colNum+1+dim(logHypChosen.GPSurvInformedARDV4.func)[2]+dim(logHypChosen.GPSurvInformedARDV4.length)[2]
	}
	if('GPSurvRF'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPNonSurvNoCens.noise
		colNum 																							<- colNum+1
	}
	if('GPSurvBIC'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvBIC.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvBIC.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvBIC.length)[2])] 			<- logHypChosen.GPSurvBIC.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvBIC.length)[2])
	}
	if('GPSurvInfMed'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvInfMed.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvInfMed.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvInfMed.length)[2])] 			<- logHypChosen.GPSurvInfMed.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvInfMed.length)[2])
	}
	if('GPSurvInfUnif'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvInfUnif.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvInfUnif.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvInfUnif.length)[2])] 		<- logHypChosen.GPSurvInfUnif.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvInfUnif.length)[2])
	}
	if('GPSurvBICForward'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvBICForward.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvBICForward.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvBICForward.length)[2])] 		<- logHypChosen.GPSurvBICForward.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvBICForward.length)[2])
	}
	if('GPSurvBICBackward'%in%models){
		logHypChosenToBoxplot[,colNum] 																	<- logHypChosen.GPSurvBICBackward.noise
		logHypChosenToBoxplot[,colNum+1] 																<- logHypChosen.GPSurvBICBackward.func
		logHypChosenToBoxplot[,(colNum+2):(colNum+1+dim(logHypChosen.GPSurvBICBackward.length)[2])] 	<- logHypChosen.GPSurvBICBackward.length
		colNum 																							<- colNum+(2+dim(logHypChosen.GPSurvBICBackward.length)[2])
	}


	# Positions to group data #
	positions 		<- numeric(length(toBoxplot))
	nextPosition 	<- 1
	for(i in 1:length(unique(toBoxplot))){
		positions[toBoxplot==unique(toBoxplot)[i]] <- nextPosition:(nextPosition+sum(toBoxplot==unique(toBoxplot)[i])-1)
		nextPosition <- nextPosition+sum(toBoxplot==unique(toBoxplot)[i])-1+3
	}

	# Colours for each model #
	# 		c('GPNonSurvNoCens','GPSurvSqExp','GPSurvARD','GPSurvInformedARD','GPSurvInformedARDV2','GPSurvInformedARDV3','GPSurvInformedARDV4','GPSurvRF','GPSurvBIC','GPSurvInfMed','GPSurvInfUnif','GPSurvBICForward','GPSurvBICBackward')
	# cols 			<- c('aquamarine3','chartreuse3','springgreen3','cadetblue3','steelblue3','cyan3','dodgerblue2','olivedrab3','turquoise3','palegreen2','seagreen3','turquoise3','turquoise4')
	cols 			<- c('chartreuse3','dodgerblue1','chartreuse3','dodgerblue1','chartreuse3','dodgerblue1','chartreuse3','dodgerblue1','chartreuse3','dodgerblue1','chartreuse3','dodgerblue1','chartreuse3')
	modelNamesPlot 	<- c('GP','GPS3SqExp','GPS3ARD','GPS3InformARDV1','GPS3InformARDV2','GPS3InformARDV3','GPS3InformARDV4','GPS3RF','GPS3BIC','GPInfMed','GPInfUnif','GPS3BICF','GPS3BICB')

	# Make boxplots #
	if(('GPSurvARD'%in%models&!separateARD)|!('GPSurvARD'%in%models)){
		if(dataSource=='Generate'){
			ylim <- c(min(c(logHypChosenToBoxplot,unlist(logHypGenerate))),max(c(logHypChosenToBoxplot,unlist(logHypGenerate))))
		} else {
			ylim <- c(min(logHypChosenToBoxplot),max(logHypChosenToBoxplot))
		}
		boxplot(logHypChosenToBoxplot,at=positions,col=cols[toBoxplot],ylim=ylim,
				main=paste0('models = ',paste0(modelNamesPlot[allModels%in%models],collapse=', ')),
				names=boxplotLabels,las=2,cex.main=0.7,cex.axis=0.7)
		mtext(text=modelNamesPlot[allModels%in%models],side=1,at=positions[!duplicated(toBoxplot)],line=3,col=cols[toBoxplot][!duplicated(cols[toBoxplot])],cex=0.7,adj=0)
		if(dataSource=='Generate'&plotHypRef){
			abline(h=logHypGenerate$noise,col='lightgrey',lty=6)
			abline(h=logHypGenerate$func,col='lightgrey',lty=4)
			abline(h=logHypGenerate$length,col='lightgrey',lty=3)
			legend('topright',legend=c(expression(paste('gen',' ',log,' ',sigma[n]^2)),expression(paste('gen',' ',log,' ',sigma[f]^2)),'log l'),col='lightgrey',pch=c(NA,NA,NA),lty=c(6,4,3),cex=0.5)
		}
	} else if('GPSurvARD'%in%models&separateARD){
		if(dataSource=='Generate'){
			ylim1 <- c(min(c(logHypChosenToBoxplot,unlist(logHypGenerate))),max(c(logHypChosenToBoxplot,unlist(logHypGenerate))))
			ylim2 <- c(min(c(cbind(logHypChosen.GPSurvARD.noise,logHypChosen.GPSurvARD.func,logHypChosen.GPSurvARD.length,logHypChosenToBoxplot),unlist(logHypGenerate))),max(c(cbind(logHypChosen.GPSurvARD.noise,logHypChosen.GPSurvARD.func,logHypChosen.GPSurvARD.length,logHypChosenToBoxplot),unlist(logHypGenerate))))
		} else {
			ylim1 <- c(min(logHypChosenToBoxplot),max(logHypChosenToBoxplot))
			ylim2 <- c(min(cbind(logHypChosen.GPSurvARD.noise,logHypChosen.GPSurvARD.func,logHypChosen.GPSurvARD.length,logHypChosenToBoxplot)),max(cbind(logHypChosen.GPSurvARD.noise,logHypChosen.GPSurvARD.func,logHypChosen.GPSurvARD.length,logHypChosenToBoxplot)))
		}
		boxplot(logHypChosenToBoxplot,at=positions,col=cols[toBoxplot],ylim=ylim1,
				main=paste0('models = ',paste0(modelNamesPlot[-3][(allModels[-3])%in%models],collapse=', ')),
				names=boxplotLabels,las=2,cex.main=0.7,cex.axis=0.7)
		mtext(text=modelNamesPlot[-3][(allModels[-3])%in%models],side=1,at=positions[!duplicated(toBoxplot)],line=3,col=cols[toBoxplot][!duplicated(cols[toBoxplot])],cex=0.7,adj=0)

		if(dataSource=='Generate'&plotHypRef){
			abline(h=logHypGenerate$noise,col='lightgrey',lty=6)
			abline(h=logHypGenerate$func,col='lightgrey',lty=4)
			abline(h=logHypGenerate$length,col='lightgrey',lty=3)
			legend('topright',legend=c(expression(paste('gen',' ',log,' ',sigma[n]^2)),expression(paste('gen',' ',log,' ',sigma[f]^2)),'log l'),col='lightgrey',pch=c(NA,NA,NA),lty=c(6,4,3),cex=0.5)
		}
		boxplot(cbind(logHypChosen.GPSurvARD.noise,logHypChosen.GPSurvARD.func,logHypChosen.GPSurvARD.length),ylim=ylim2,
				col=cols[3],main='model = GPSurvARD',names=c('noise','func',paste0('length',1:dim(logHypChosen.GPSurvARD.length)[2])),las=2,
				cex.main=0.7,cex.axis=0.7,medlwd=0.5,boxlwd=0.5,whisklwd=0.5,staplelwd=0.5)
		if(dataSource=='Generate'&plotHypRef){
			abline(h=logHypGenerate$noise,col='lightgrey',lty=6)
			abline(h=logHypGenerate$func,col='lightgrey',lty=4)
			abline(h=logHypGenerate$length,col='lightgrey',lty=3)
			legend('topright',legend=c(expression(paste('gen',' ',log,' ',sigma[n]^2)),expression(paste('gen',' ',log,' ',sigma[f]^2)),'log l'),col='lightgrey',pch=c(NA,NA,NA),lty=c(6,4,3),cex=0.5)
		}
	}

}