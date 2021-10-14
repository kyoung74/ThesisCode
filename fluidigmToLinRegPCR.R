fluidigmToLinRegPCR=function(rawSYBRGreen,bkgdSYBRGreen,rawROX,bkgdROX,samples,assays){
	# Input files are .csv exported from Fluidigm software from .bml to Table with Raw Data. Each .csv contains the data from all sample-assay combinations. samples is a .csv of sample IDs under heading "SampleID" (include Pos and Neg controls) and assays is a .csv of primer probes under heading "Primer".
	library(stringr)
	
	# Read .csv files
	rawSYBRGreen=read.csv(rawSYBRGreen)
	bkgdSYBRGreen=read.csv(bkgdSYBRGreen) 
	rawROX=read.csv(rawROX)
	bkgdROX=read.csv(bkgdROX)
	samples=read.csv(samples)
	assays=read.csv(assays)
	
	# Calculate Rn
	Rn=(rawSYBRGreen[,-1]-bkgdSYBRGreen[,-1])/(rawROX[,-1]-bkgdROX[,-1])
	
	# Match samples and assays to rows
	cycleNum=length(Rn)
	Rn$ChamberID=rawSYBRGreen[,1]
	names(Rn)=c(1:cycleNum,"ChamberID")
	saNames=str_match(Rn$ChamberID,"(S)(.{2})(-A)(.{2})$")
	Rn$S=as.numeric(as.character(saNames[,3]))
	Rn$A=as.numeric(as.character(saNames[,5]))
	Rn$S=factor(Rn$S)
	Rn$A=factor(Rn$A)
	levels(Rn$S)=samples$SampleID
	levels(Rn$A)=assays$Primer
	
	# Creating LinRegPCR Excel Sheet
	data=c()
	dataFrame=c()
	numChambers=dim(Rn)[1]
	for(n in 1:numChambers){
		print(paste(n," of ",numChambers,sep=""))
		data$Well=rep(Rn$S[n],cycleNum)
		data$Cycle=c(1:cycleNum)
		data$TargetName=rep(Rn$A[n],cycleNum)
		RnFrame=data.frame(t(Rn[n,1:cycleNum]))
		data$Rn=RnFrame[,1]
		dataFrame=rbind(dataFrame,data.frame(data))
	}
	
	# Sorting data alphabetically and Splitting into manageable subgroups
	sortedData=dataFrame[order(dataFrame$TargetName),]
	dataLength=dim(sortedData)[1]
	sampleNum=dim(samples)[1]
	twelfth=dataLength/(cycleNum*sampleNum*12)
	end=cycleNum*sampleNum*12
	finalFrame=sortedData[1:(cycleNum*sampleNum*12),]
	start=(cycleNum*sampleNum*12)+1
	for(n in 2:twelfth){
		end=cycleNum*sampleNum*12*n
		finalFrame=cbind(finalFrame,sortedData[start:end,])
		start=end+1
	}
	write.csv(finalFrame,"FileForLinRegPCR.csv")
}