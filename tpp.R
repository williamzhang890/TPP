init = function(frame, identifierColumn, coverageColumn, numPeptidesColumn, qpsmColumn,
	vehicle1Columns, treatment1Columns, vehicle2Columns = NULL, treatment2Columns = NULL, 
	experimentNames = c("Vehicle1", "Treatment1", "Vehicle2", "Treatment2"), 
	temperatures, tagNames, resultPath = getwd()) {
	
	# checks if the TPP package has been loaded in and loads it in if it hasn't been loaded
	if (!require("TPP",character.only = TRUE)) {
		source("https://bioconductor.org/biocLite.R");
		biocLite("TPP");
		library("TPP");
	}

	# checks if directory where files will be saved exists
	if (!dir.exists(resultPath)) {	
		dir.create(resultPath);
	}

	# function used to convert quantifications to a relative format
	toRelative = function(columns) {
		return (function(row) {
			row = as.numeric(row);
			row = row[columns] / row[columns[1]];
			return (row);
		});
	}

	# generates qssm values
	coverage = as.numeric(frame[,coverageColumn]);
	peptides = as.numeric(frame[,numPeptidesColumn]);
	qpsm = as.numeric(frame[,qpsmColumn]);
	qupm = peptides * coverage/100;
	qssm = qpsm / qupm;
	frame[, "qssm"] = qssm;
	qssmIndex = which(colnames(frame) == "qssm");

	# column names used by TPP
	colNames = c("Accession", "qssm", lapply(tagNames, 
		function(s) {
			return (paste("rel_fc_",s, sep = ""));			
		}));

	# creates config file and creates data frames that will be used by the TPP package
	config = data.frame(matrix(vector(), 0, 14), stringsAsFactors = F, check.names = F);
	names(config) = c("Experiment", "Condition", "ComparisonVT1", "ComparisonVT2", tagNames);

	# vehicle 1
	config[1,] = c(experimentNames[1], "Vehicle", "X", "", temperatures); # row 1 of config
	frame[, vehicle1Columns] = suppressWarnings(t(apply(frame, 1, toRelative(vehicle1Columns)))); # converts all vehicle 1 columns to a normalized format 
	vehicle1 = frame[, c(identifierColumn, qssmIndex, vehicle1Columns)]; # creates data frame for vehicle 1
	names(vehicle1) = colNames;
	write.csv(vehicle1, file = file.path(resultPath,"vehicle1.csv"), row.names=F);
	
	# treatment 1
	config[2,] = c(experimentNames[2], "Treatment", "X", "", temperatures);	# row 2 of config
	frame[, treatment1Columns] = suppressWarnings(t(apply(frame, 1, toRelative(treatment1Columns)))); # converts all treatment 1 columns to a normalized format 
	treatment1 = frame[, c(identifierColumn, qssmIndex, treatment1Columns)]; # creates data frame for treatment 1
	names(treatment1) = colNames;
	write.csv(treatment1, file = file.path(resultPath,"treatment1.csv"), row.names=F);


	if (!is.null(vehicle2Columns) && !is.null(treatment2Columns)) { # checks if a parallel trials were run
		# vehicle 1
		config[3,] = c(experimentNames[3], "Vehicle", "X", "", temperatures); # row 3 of config
		frame[, vehicle2Columns] = suppressWarnings(t(apply(cetsa, 1, toRelative(vehicle2Columns)))); # converts all vehicle 2 columns to a normalized format 
		vehicle2 = cetsa[, c(identifierColumn, qssmIndex, vehicle2Columns)]; # creates data frame for vehicle 2
		names(vehicle2) = colNames;
		write.csv(vehicle2, file = file.path(resultPath,"vehicle2.csv"), row.names=F);

		# treatment 1
		config[4,] = c(experimentNames[4], "Treatment", "X", "", temperatures);	# row 4 of config
		frame[, treatment2Columns] = suppressWarnings(t(apply(frame, 1, toRelative(treatment2Columns)))); # converts all treatment 2 columns to a normalized format 
		treatment2 = frame[, c(identifierColumn, qssmIndex, treatment2Columns)]; # creates data frame for treatment 2
		names(treatment2) = colNames;
		write.csv(treatment2, file = file.path(resultPath,"treatment2.csv"), row.names=F);
	}
	write.csv(config, file = "config.csv", row.names = F);
}

analyze = function(configPath = "config.csv", vehicle1Path = "vehicle1.csv", treatment1Path = "treatment1.csv", 
	vehicle2Path = NULL, treatment2Path = NULL, experimentNames = c("Vehicle1", "Treatment1", "Vehicle2", "Treatment2"), 
	resultPath = getwd()) {
	
	# checks if the TPP package has been loaded in and loads it in if it hasn't been loaded
	if (!require("TPP",character.only = TRUE)) {
		source("https://bioconductor.org/biocLite.R");
		biocLite("TPP");
		library("TPP");
	}

	# checks if directory where files will be saved exists
	if (!dir.exists(resultPath)) {	
		dir.create(resultPath);
	}
	
	# load in config, vehicle1, and vehicle2 dataframes
	config = read.csv(configPath, check.names = F);
	vehicle1 = read.csv(vehicle1Path);
	treatment1 = read.csv(treatment1Path);

	# instantiate parameters used by TPP
	data = list(vehicle1, treatment1);
	names(data) = experimentNames[1:2];
	
	# load in vehicle2 and treatment2 dataframes and add them to data
	if (!is.null(vehicle2Path) && !is.null(treatment2Path)) {
		vehicle2 = read.csv(vehicle2Path);
		treatment2 = read.csv(treatment2Path);
		data = c(data, list(vehicle2, treatment2));
		names(data)[3:4] = experimentNames[3:4];
	}

	# runs TPP's analysis but does not plot the curves
	results = analyzeTPPTR(configTable = config, methods = "meltcurvefit", data, resultPath = resultPath, plotCurves = F, idVar = "Accession");
	write.csv(results, "results.csv", row.names=F);
}

getSignificant = function(resultPath = "results.csv", experimentNames = c("Vehicle1", "Treatment1", "Vehicle2", "Treatment2")) {
	
	results = read.csv(resultPath, check.names = F);

	# load in comparison 1 indexes
	pVal1 = which(names(results) == paste("pVal_adj_", experimentNames[2], "_vs_", experimentNames[1], sep = ""));
	minSlope1 = which(names(results) == paste("min_slope_", experimentNames[2], "_vs_", experimentNames[1], sep = ""));
	r2Vehicle1 = which(names(results) == paste("R_sq_", experimentNames[1], sep = ""));
	r2Treatment1 = which(names(results) == paste("R_sq_", experimentNames[2], sep = ""));
	plateau1 = which(names(results) == paste("plateau_", experimentNames[1], sep = ""));
	
	fulfillsAll = which(names(results) == "fullfills_all_4_requirements");
	
	# apply filters
	if (length(fulfillsAll) == 0) { # check if parallel trials haven't been ran
		# p-value < .05
		# pg 9 of Savitski et al.
		results = results[which(results[,pVal1] < .05),];
		
		# min slope < -.06 
		# pg 7 of TPP documentation
		results = results[which(results[,minSlope1] < -.06),]; 
		
		# r-squared > .8 
		# pg 9 of Savitski et al.
		results = results[which(results[,r2Vehicle1] > .8),];
		results = results[which(results[,r2Treatment1] > .8),];

		# plateau for control < .3
		# pg 9 of Savitski et al.
		results = results[which(results[,plateau1] < .3),]; 
	} else {
		# load in comparison 2 indexes
		r2Vehicle2 = which(names(results) == paste("R_sq_", experimentNames[3], sep = ""));
		r2Treatment2 = which(names(results) == paste("R_sq_", experimentNames[4], sep = ""));
		plateau2 = which(names(results) == paste("plateau_", experimentNames[3], sep = ""));
		
		# min slope < -.06
		# melting point differences in the vehicle vs treatment > melting point difference between the two vehicles
		# one of the p values < 0.05 and the other one < 0.1
		# melting point shifts have the same sign
		# pg 7 of TPP documentation
		results = results[which(fulfillsAll),]; 
		
		# one r-squared > .9 and the other r-squared > .8
		# pg 9 of Savitski et al.
		results = results[which(results[,r2Vehicle1] > .8 & results[,r2Vehicle2] > .8),];
		results = results[which(min(results[,r2Vehicle1], results[,r2Vehicle2]) > .9),];
		results = results[which(results[,r2Treatment1] > .8 & results[,r2Treatment2] > .8),];
		results = results[which(min(results[,r2Treatment1], results[,r2Treatment2]) > .9),];

		# plateaus for controls < .3
		# pg 9 of Savitski et al.
		results = results[which(results[,plateau1] < .3 & results[,plateau2] < .3)]
	}

	return (results$Protein_ID);
}

plotCurves = function(configPath = "config.csv", vehicle1Path = "vehicle1.csv", treatment1Path = "treatment1.csv", 
	vehicle2Path = NULL, treatment2Path = NULL, experimentNames = c("Vehicle1", "Treatment1", "Vehicle2", "Treatment2"), 
	resultPath = getwd(), proteinList) {

	# checks if the TPP package has been loaded in and loads it in if it hasn't been loaded
	if (!require("TPP",character.only = TRUE)) {
		source("https://bioconductor.org/biocLite.R");
		biocLite("TPP");
		library("TPP");
	}

	# checks if directory where files will be saved exists
	if (!dir.exists(resultPath)) {	
		dir.create(resultPath);
	}

	# load in config, vehicle1, and vehicle2 dataframes
	config = read.csv(configPath, check.names = F);
	vehicle1 = read.csv(vehicle1Path);
	treatment1 = read.csv(treatment1Path);

	# instantiate parameters used by TPP
	data = list(vehicle1, treatment1);
	names(data)[1:2] = experimentNames[1:2];
	
	# load in vehicle2 and treatment2 dataframes and add them to data
	if (!is.null(vehicle2Path) && !is.null(treatment2Path)) {
		vehicle2 = read.csv(vehicle2Path);
		treatment2 = read.csv(treatment2Path);
		data = c(data, list(vehicle2, treatment2));
		names(data)[3:4] = experimentNames[3:4];
	}
	
	# normalizes data
	dataNormalized = tpptrNormalize(data=tpptrImport(configTable = config, data = data, idVar = "Accession"))[["normData"]];

	print(dataNormalized)
	# plot curves for proteins in proteinList
	tpptrCurveFit(data = lapply(dataNormalized, function(d)
		d[Biobase::featureNames(d) %in% proteinList,]), resultPath = resultPath);
}

poissonFit = function(frame, identifierColumn, vehicleColumns, treatmentColumns, temperatures, resultPath = getwd()) {
	
	# essential functions
	poisson = function(x,L,a,h,C) { # function used to model the poisson curve
		x[which(x<h)] = h;
		return(((1-C) * L ^ (a * (x - h)))/gamma(a * (x - h) + 1) + C);
	}
	angle = function(a, b) { # function used to calculate the angle between 2 vectors
		# cos(theta) = (a dot b) / (||a|| * ||b||)
		p = a * b;
		dot = sum(na.omit(p));
		v1 = a[which(!is.na(p))];
		v2 = b[which(!is.na(p))];
		v1 = sqrt(sum(v1 * v1));
		v2 = sqrt(sum(v2 * v2));
		return (acos(dot/(v1*v2)));
	}
	model = function(x,y) { # function used to generate the model
		m = try(nls(y ~ poisson(x,L,a,h,C),
			start = c(
				L = 1.5,
				a = 1/5,
				h = 40,
				C = 0),
			algorithm = "port",
			lower = c(
				L = 0,
				a = 0,
				h = 37,
				C = 0),
			upper = c(
				L = 10,
				a = 1,
					h = 70,
				C = .5)),silent = T); # try suppresses warning if nls failed
		return (m);
	}
	deriv = function(func, x) { # calculates the derivative with respect to x of a function at a point
		h = .001;
		return((func(x + h) - func(x))/h);
	}

	# essential columns
	names = frame[,identifierColumn];
	vehicle = frame[,vehicleColumns];
	treatment = frame[,treatmentColumns];
	
	folderPath = file.path(resultPath, "Melting_Curves");
	# checks if directory where files will be saved exists
	if (!dir.exists(folderPath)) {	
		dir.create(folderPath);
	}
	
	# columns of output table

	angularSimilarity = numeric(nrow(frame));
	
	# parameters for vehicle
	vehicle_Fitted = logical(nrow(frame));
	vehicle_L = numeric(nrow(frame));
	vehicle_a = numeric(nrow(frame));
	vehicle_h = numeric(nrow(frame));
	vehicle_C = numeric(nrow(frame));
	vehicle_slopeAtMid = numeric(nrow(frame));
	vehicle_meltingPoint = numeric(nrow(frame));
	
	# parameters for treatment
	treatment_Fitted = logical(nrow(frame));
	treatment_L = numeric(nrow(frame));
	treatment_a = numeric(nrow(frame));
	treatment_h = numeric(nrow(frame));
	treatment_C = numeric(nrow(frame));
	treatment_slopeAtMid = numeric(nrow(frame));
	treatment_meltingPoint = numeric(nrow(frame));

	meltingShifts = numeric(nrow(frame));

	for (i in 1:nrow(frame)) {
		# converts quantifications to relative format
		v = as.numeric(vehicle[i,] / vehicle[i,1]);
		t = as.numeric(treatment[i,] / treatment[i,1]);
		angularSimilarity[i] = 1 - angle(v, t)/pi;

		# temperature vectors that do not contain temperatures that correspond to missing values
		vt = temperatures[which(!is.na(v))];
		tt = temperatures[which(!is.na(t))];
		# delete missing values from quantifications
		v = na.omit(v);
		t = na.omit(t);

		# set default for every column to NA
		vehicle_L[i] = NA;
		vehicle_a[i] = NA;
		vehicle_h[i] = NA;
		vehicle_C[i] = NA;
		vehicle_meltingPoint[i] = NA;
		vehicle_slopeAtMid[i] = NA;
		treatment_L[i] = NA;
		treatment_a[i] = NA;
		treatment_h[i] = NA;
		treatment_C[i] = NA;
		treatment_meltingPoint[i] = NA;
		treatment_slopeAtMid[i] = NA;
		meltingShifts[i] = NA;

		if (length(v) > 0 && length(t) > 0) { # checks if a point exists so that a plot can be made
			# generate poisson models
			vm = model(vt, v);
			tm = model(tt, t);

			# plot scatterplot
			pdf(file.path(folderPath,paste(names[i],".pdf", sep = "")));
			plot(vt, v, ylim = c(0, max(v,t)),
				main = names[i], ylab = "Relative Protein Amount", xlab = "Temperature (Celsius)");
			points(tt, t, pch = 2, col = "green");
			if(is.list(vm)) { # checks if poisson model was fitted for the vehicle
				vehicle_Fitted[i] = T;
				
				# save graph parameters
				cs = as.numeric(coef(vm));
				vehicle_L[i] = cs[1];
				vehicle_a[i] = cs[2];
				vehicle_h[i] = cs[3];
				vehicle_C[i] = cs[4];

				# plot the model's curve
				fun = function(x){return (poisson(x,cs[1],cs[2],cs[3],cs[4]))}; # function used to plot the curve
				plot(fun,
					from = 37, to = 70, lty = "dashed", add = T);

				# find melting point
				mp = as.numeric(optimize(function(x){return (abs(fun(x) - .5))}, interval = c(37,70)))[1];
				if (abs(fun(mp)-.5) < .01) { # checks if mp is really the melting point and not just one of the endpoints
					vehicle_meltingPoint[i] = mp;
					vehicle_slopeAtMid[i] = deriv(fun, mp);
					points(mp, fun(mp), pch = 4, cex = 2); # draws X at melting point
				}
			} 

			if(is.list(tm)) { # checks if poisson model was fitted for the treatment
				treatment_Fitted[i] = T;

				# save graph parameters
				cs = as.numeric(coef(tm));
				treatment_L[i] = cs[1];
				treatment_a[i] = cs[2];
				treatment_h[i] = cs[3];
				treatment_C[i] = cs[4];

				# plot the model's curve
				fun = function(x){return (poisson(x,cs[1],cs[2],cs[3],cs[4]))};
				plot(fun,
					from = 37, to = 70, lty = "dashed", col = "green", add = T);

				# find melting point
				mp = as.numeric(optimize(function(x){return (abs(fun(x) - .5))}, interval = c(37,70)))[1];
				if (abs(fun(mp)-.5) < .01) { # checks if mp is really the melting point and not just one of the endpoints
					treatment_meltingPoint[i] = mp;
					treatment_slopeAtMid[i] = deriv(fun, mp);
					points(mp, fun(mp), pch = 4, col = "green", cex = 2); # draws X at melting point
				}
			}

			meltingShifts[i] = treatment_meltingPoint[i] - vehicle_meltingPoint[i]; # calculates melting point shifts. if one of the points was NA, it will still be NA
			legend("topright", pch = c(1, 2), lty = "dashed", col = c("black", "green"), legend = c("Control", "Treated"));
			dev.off(); # save pdf
		}
	}

	write.csv(data.frame(names, angularSimilarity, vehicle_Fitted, vehicle_L, vehicle_a, vehicle_h, vehicle_C, vehicle_meltingPoint, vehicle_slopeAtMid,
						treatment_Fitted, treatment_L, treatment_a, treatment_h, treatment_C, treatment_meltingPoint, treatment_slopeAtMid,
				   meltingShifts),
		file.path(resultPath, "results.csv"), row.names = F); # save results
}