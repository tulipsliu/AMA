callAMA <-
function(cofh, neq, leads, lags, qmax = -42){ #qmax is initialized to 150000.  To change, use "qmax = <value>" as the last argument.
	#make sure inputs are okay.  initialize res to "Error".
	res = "Error: could not compute"
	hrows = neq;
	hcols = neq*(leads + lags + 1)
	cofh = as.matrix(cofh)
	cofh = matrix(cofh, hrows)
	if(qmax == -42){
		qmax = (neq * (lags + leads))*(neq * (lags + leads))
	}

	if(nrow(cofh) != hrows){
		print("Error: neq not equal to number of rows actually in H matrix.")
	}

	if(ncol(cofh) != hcols){
		print("Error: neq*(leads + lags + 1) not equal to number of columns actually in H matrix.")
	}
	
	else if(hrows != as.integer(hrows) || hrows <= 0){
		print("Error: hrows must be a positive integer.")
	}

	else if(hcols != as.integer(hcols) || hrows <= 0){
		print("Error: hcols must be a positive integer.")	
	}

	else if(neq != as.integer(neq) || neq <= 0){
		print("Error: neq must be a positive integer.")
	}

	else if(leads != as.integer(leads) || leads <= 0){
		print("Error: leads must be a positive integer.")
	}

	else if(lags != as.integer(lags) || lags <= 0){
		print("Error: lags must be a positive integer.")
	}

	else if(qmax != as.integer(qmax) || qmax <= 0){
		print("Error: qmax must be a positive integer.")
	}

	#and if they're okay, call AMA (In the C lingo, this is called SparseAIM)
	else{

		Rcofh = as.double(cofh)
		Rhrows = as.integer(hrows)
		Rhcols = as.integer(hcols)
		Rneq = as.integer(neq)
		Rleads = as.integer(leads)
		Rlags = as.integer(lags)
		Rnstate = as.integer(leads + lags + 1)
		Rqmax = as.integer(qmax)
		Rptr = as.double(0)
		Rcofb = as.double(mat.or.vec((neq*leads) * (neq*lags), 1))
		RcofQ = as.double(mat.or.vec((neq*leads) * (neq*(lags+leads)), 1))

		res <- .C("callSparseAimFromR", Rcofh, Rhrows, Rhcols, Rneq, 
			Rleads, Rlags, Rnstate, Rqmax, Rptr, Rcofb, RcofQ)
	}

	res

} #end of main wrapper function



#---------------GENERATES B (REDUCED FORM COEFFICIENTS) MATRIX------------------------#

genBmat <- function(cofh, neq, leads, lags, qmax = -42){
	if(qmax == -42){
		qmax = (neq * (lags + leads))*(neq * (lags + leads))
	}
	hrows = neq;
	hcols = neq*(leads + lags + 1)
	res <- callAMA(cofh, neq, leads, lags)
	B = res[10][[1]]
	B = as.matrix(B)
	B = matrix(B, nrow=neq,ncol=(neq*lags),byrow=FALSE )
	B

}

#-------------------------------------END genBmat--------------------------------------#


#---------------GENERATES Q (ASYMPTOTIC LINEAR CONSTRAINTS) MATRIX----------------------#

genQmat <- function(cofh, neq, leads, lags, qmax = -42){
	if(qmax == -42){
		qmax = (neq * (lags + leads))*(neq * (lags + leads))
	}
	hrows = neq;
	hcols = neq*(leads + lags + 1)
	res <- callAMA(cofh, neq, leads, lags)
	Q = res[11][[1]]
	Q = as.matrix(Q)
	Q = matrix(Q, (neq*leads) )
	Q
}
#-------------------------------------END genQmat--------------------------------------#





#-----------shift right function--used in genScof------------#
SPshiftright <- function(x, n){
	rows <- nrow(x);
	cols <- ncol(x);
	z = mat.or.vec(rows, n);
	y = cbind(z,x);
	y = y[,1:cols]
	y;
}
#-------------------------------------END SPshiftright---------------------------------#




#--------------GENERATES S (OBSERVABLE STRUCTURE) MATRIX------------------#
genScof <- function(hmat, bmat, neq, leads, lags){

	#--------initialize-----------
	nlag = lags;
	nlead = leads;
	hmat2 = as.matrix(hmat);
	hmat2 = matrix(hmat2, neq)
	bmat2 = cbind(bmat, -diag(neq))
	scof = mat.or.vec(neq, neq*(nlag+1))
	q = mat.or.vec(neq*nlead, neq*(nlag+nlead))
	rc = nrow(bmat2)
	cc = ncol(bmat2)
	q[1:rc, 1:cc] = bmat2;
	qcols = neq*(nlag + nlead)

	if(nlead > 1){ #shift right
		for(i in 1:(nlead-1)){
			rows = i*neq + (1:neq)
			q[rows, ] = SPshiftright(q[rows-neq,], neq) 
		}
	} #end if-statement

	L = (1:(neq*nlag))
	R = ((neq*nlag+1) : (neq*(nlag+nlead)))
	
	#----------here we solve the linear system.-------------
	qRight = (as.matrix(-q[,R]));
	qLeft = (as.matrix(q[,L]));

	q[,L] = solve(qRight, qLeft);

	#call sparskit
	#qRightSparse <- as.matrix.csr(qRight)
	#qLeftSparse <- as.matrix.csr(qLeft)

	#----------put it all together---------------
	minus = 1 : (neq*(nlag+1));
	plus = (neq*(nlag+1)+1) : (neq*(nlag+1+nlead));
	scof[,(neq+1):(neq*(nlag+1))] = hmat2[,plus]%*%q[,L]
	scof = scof + hmat2[,minus] 
	scof	

}
#-----------------------------END genScof------------------------------------


#-------------------------GET RETURN CODE FOR ERRORS-------------------------

getReturnCode <- function(sparseAimObject){
	returnInt = sparseAimObject[9][[1]]
	sprintf("Return Code: %d", returnInt)

	if(returnInt == 0){
		returnCode = "AIM gives unique solution."
	}

	else if(returnInt == 2000){
		returnCode = "Stacked system of H matrices not full rank"
	}
	
	else if(returnInt == 2001){	 
		returnCode = "sparseAim_PRECONDITIONS_VIOLATED: Inappropriate dimensions for H/Q matrices."
	}		

	else if(returnInt == 2002){	
		returnCode = "autoRegression_POSTCONDITIONS_VIOLATED: unconstrained autoregression returns invalid matrices"
	}		

	else if(returnInt == 2003){	
		returnCode = "augmentQmatWithInvariantSpaceVectors_PRECONDITIONS_VIOLATED: need positive number of constraints and nonnegative number of auxiliary initial conditions"
	}		


	else if(returnInt == 2004){	
		returnCode = "augmentQmatWithInvariantSpaceVectors_POSTCONDITIONS_VIOLATED: Q matrix has invalid dimension."
	}		

	else if(returnInt == 2005){	#
		returnCode = "shiftRightAndRecord_PRECONDITIONS_VIOLATED; too many exact or numeric shiftrights."
	}		

	else if(returnInt == 2006){	
		returnCode = "annihilateRows_POSTCONDITIONS_VIOLATED: row annihilator yields invalid matrix"
	}		

	else if(returnInt == 2007){	
		returnCode = "Calculations require larger number for maxNumberOfHElements (qmax)"
	}		

	else if(returnInt == 2008){	
		returnCode = "Matrix A too large"
	}	

	else if(returnInt == 2009){	
		returnCode = "Not enough large roots"
	}	

	else if(returnInt == 2010){	
		returnCode = "Too many large roots"
	}		

	else{
		returnCode = "Error: return code not properly specified."
	}

	returnCode
}
#---------------------------END getReturnCode--------------------------------



#---------------------------COMPUTE STOCHASTIC TRANSITION MATRICES-----------

getStochTrans <- function(hmat, scof){ #hmat and scof must be in matrix format, i.e. have appropriate # rows & cols
	
	hrows = nrow(hmat);
	hcols = ncol(hmat);
	srows = nrow(scof);
	scols = ncol(scof);
	neq = hrows;
	nlags = (scols/neq) - 1	
	sPartial = scof[,(nlags*neq + 1):((nlags+1)*neq)]
	sinv = solve(sPartial)  #invert the right block of the matrix S
	if(nlags > 1){
		zer0s = mat.or.vec( ( (nlags-1)*neq ),neq )
		iden = diag((nlags-1)*neq)
		sRow = -sinv %*% scof[,1:(nlags*neq)]
		aTrans = cbind(zer0s, iden);
		aTrans = rbind(aTrans, sRow);
		bTrans = cbind(zer0s, sinv);
	}
	else{
		aTrans = -sinv %*% scof[,(1:neq)]
		bTrans = sinv;
	}

	returnObject = list(aTrans, bTrans);
	returnObject
}

#-------------------------END getStochTrans---------------------------------


#------------------COMPUTE INHOMOGENEOUS MATRICES (phi, F, theta)--------

getFactorMatrices <- function(hmat, bmat, neq, leads, lags){ 
	hmat0 = hmat[, (lags*neq+1):((lags+1)*neq)]
	hcols = ncol(hmat)
	bcols = ncol(bmat)
	bmatRight = bmat[,bcols-neq + 1:bcols]

	#------compute phi--------------
	hmatplus = hmat[,((lags+1)*neq + 1):hcols ]
	phiInv = hmat0 + hmatplus %*% bmatRight
	phi = solve(phiInv)
	
	#-------compute F---------------
	phihp = -phi%*% hmatplus
	if(leads > 1){
		topMatrix = mat.or.vec(neq*(leads-1), neq*leads)
		bottomMatrix = NULL;
		for(i in 1:leads){
			j = leads - i;
			zeroAndIden = mat.or.vec(neq*j, neq)
			zeroAndIden = rbind(zeroAndIden, diag(neq))
			bPortion = NULL
			if(i != 1){
				bPortion = bmatRight[1:((i-1)*neq ),]
				bPortionSparse = phihp %*% bPortion
				topMatrix[, ((i-1)*neq + 1):(i*neq)] <- diag(neq)

			}
			nextMat = rbind(zeroAndIden, bPortion)
			cbind(bottomMatrix, nextMat)
		} #end for-loop

		fmat = rbind(topMatrix, bottomMatrix)
	} #end if-statement	

	else{
		fmat = phihp		
	}

	outList <- list(phi, fmat)
	outList

}
#-----------------------END getFactorMatrices-------------------------------------




#----------------------PARSE INPUT FROM MODELEZ SYNTAX----------------------------

#returns filenames of R scripts which can then be called to load the matrices
genHmatrixScripts <-function(modelFileNameFull){ 
#	javaIOstream <- new(J("java.io.FileInputStream"), modelFileName)
#	modelFileName = strsplit(modelFileNameFull, "")[[1]][1]
	modelFileName = modelFileNameFull
	hungryCaterpillar <- new(J("modelez_R.RunTheParser"))
	javaFileName = .jnew("java/lang/String", modelFileName);	#create Java string object to pass to runParser
#	fileStringArray = .jcast(javaFileName, "[Ljava/lang/String;") #and cast it as a String[]-like object
	.jcall(hungryCaterpillar, returnSig = "S", method = "runTheParser", javaFileName) 
	dataFileName = paste(modelFileName, "_SparseAimData.r", sep = "")
	matricesFileName = paste(modelFileName, "_SparseAimMatrices.r", sep = "")
#	matricesCode <- read.fwf("matricesFileName", widths = 300)
	returnObject <- list(dataFileName, matricesFileName)
	returnObject;

}

#----------------------END genHmatrixScripts------------------------------------------



#-------------PARSE INPUT FROM MODELEZ SYNTAX AND GENERATE H MATRICES-----------------

genHmat <- function(modelFileNameFull, params, wantParamVec = FALSE){
	genHmatrixScripts(modelFileNameFull)
#	modelFileName = strsplit(modelFileNameFull, "")[[1]][1]
	modelFileName = modelFileNameFull
	if(class(params) == "character"){ 
		paramTable <- read.fwf(params, widths = 250)
		rows = nrow(paramTable)
		paramVec = mat.or.vec(rows,1)
		for(i in 1:rows){
			str2 = paramTable[[1]][i];
			eval(parse(text = str2))
			txt = as.character(paramTable[[1]][i])
			splitted <- strsplit(txt, "=")
			compatibleNumber <- as.numeric(strsplit(splitted[[1]][2], ";")[1])
			paramVec[i] <- compatibleNumber
			names(paramVec)[i] <- splitted[[1]][1]
		}
	}

	else{	#i.e. if class(params) is either matrix or numeric 
		nParams <- length(params)
		for (i in 1:nParams){
			asStr = paste(names(params[i]), "=", as.character(params[i]))
			eval(parse(text = asStr))
		}
	}
cofh=0#cofh will actually be defined by sourceFileName 
#this line eliminates undefined complaint from R CMD check
	sourceFileName = paste(basename(modelFileName), "_SparseAimMatricesR.r", sep = "")
	sourceFileTable <- read.fwf(sourceFileName, widths = 250)
	rowTab = nrow(sourceFileTable)
	for(j in 1:rowTab){
		str3 = sourceFileTable[[1]][j]
		eval(parse(text = str3))
	}

	if(wantParamVec == FALSE){
		return(cofh); #sourceFile generates the H matrix
	}
	else{
		outlist <- list(cofh, paramVec)
		return(outlist)
	}
}


