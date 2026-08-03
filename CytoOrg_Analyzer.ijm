// =======================================================
// CytoOrg Analyzer
// Automated cytoskeletal organization analysis pipeline
// Fiji/ImageJ macro
//
// Version: 1.0
// Author: Perina Edoardo
// Date: 2026
// =======================================================

// ==============================
// CLICK RUN TO START
// ==============================

dir = getDirectory("Choose input folder");

// ==============================
// INSERT PARAMETERS
// ==============================

micro = fromCharCode(181);
sup2 = fromCharCode(178);

Dialog.create("Cytoskeletal Organization Analysis");

Dialog.addNumber("Image pixel size (" + micro + "/pixel)", 0.0547606);

Dialog.addCheckbox("Measure filament thickness", true);
Dialog.addSlider("Minimum intensity (0 = Auto)", 0, 255, 40);
Dialog.addNumber("Ignore thickness values below (" + micro + "m)", 0.2);
Dialog.addMessage("Thickness map display range (" + micro + "m)");
Dialog.addNumber("Minimum", 0);
Dialog.addNumber("Maximum", 1.5);
Dialog.addMessage(" ");

Dialog.addNumber("Tubeness filter sigma", 0.3);
Dialog.addCheckbox("Estimate sigma from filament thickness", false);
Dialog.addMessage(" ");

Dialog.addNumber("Exclude filaments shorter than (" + micro + "m)", 0.2);
Dialog.addNumber("Short / Intermediate cutoff (" + micro + "m)", 1.5);
Dialog.addNumber("Intermediate / Long cutoff (" + micro + "m)", 3.0);
Dialog.addCheckbox("Export measurements for individual filaments", false);
Dialog.addCheckbox("Export measurements for grouped filament", false);
Dialog.addCheckbox("Save Skeleton improved for visualization", false);
Dialog.addCheckbox("Save Branching map", false);
Dialog.addCheckbox("Save Branching map improved for visualization", false);
Dialog.addMessage(" ");

Dialog.addNumber("OrientationJ window sigma", 2);
methods = newArray(
    "Cubic Spline",
    "Finite Difference",
    "Fourier",
    "Riesz Filters",
    "Gaussian",
    "Hessian"
);
Dialog.addChoice("Gradient calculation method", methods, "Cubic Spline");
Dialog.addSlider("Direction threshold (degrees)", 1, 90, 45);
Dialog.addCheckbox("Save Orientation map", false);
Dialog.addCheckbox("Save Orientation map improved for visualization", false);
Dialog.addMessage(" ");

Dialog.addString("Choose output folder name", "Results");

Dialog.show();


// ==============================
// READ PARAMETERS
// ==============================

pixelSize = Dialog.getNumber();

measureThickness = Dialog.getCheckbox();
IntensityThreshold = round(Dialog.getNumber());
minBranchThickness = Dialog.getNumber();
thicknessMin = Dialog.getNumber();
thicknessMax = Dialog.getNumber();

sigmaTubeness = Dialog.getNumber();
automaticSigma = Dialog.getCheckbox();

minBranchLength = Dialog.getNumber();
ShortIntermediate = Dialog.getNumber();
IntermediateLong = Dialog.getNumber();
saveIndividualFilaments = Dialog.getCheckbox();
saveGroupsFilaments = Dialog.getCheckbox();
saveEnhancedSkeleton = Dialog.getCheckbox();

saveBranchMap = Dialog.getCheckbox();
saveEnhancedBranchMap = Dialog.getCheckbox();

orientationSigma = Dialog.getNumber();
gradientMethod = Dialog.getChoice();

directionChangeThreshold = Dialog.getNumber();

saveOrientationMap = Dialog.getCheckbox();
saveEnhancedOrientationMap = Dialog.getCheckbox();

OutputFolderName = Dialog.getString();

// ======================================
// GET LIST OF IMAGES IN INPUT FOLDER
// ======================================

files = getFileList(dir);
validFiles = newArray(0);

for (i = 0; i < files.length; i++) {
    name = files[i];
    if (endsWith(name, ".tif") || 
        endsWith(name, ".tiff") ||
        endsWith(name, ".png") ||
        endsWith(name, ".jpg") ||
        endsWith(name, ".jpeg")) {
        validFiles = Array.concat(validFiles, name);
    }
}

if (validFiles.length == 0) {
    exit("No valid image files found in selected folder.");
}

// ======================================
// CREATE OUTPUT FOLDER
// ======================================

if (OutputFolderName=="") {
    OutputFolderName = "Results";
}

outputDir = dir + OutputFolderName + File.separator;

if (!File.exists(outputDir)) {
    File.makeDirectory(outputDir);
}


// Hide image windows and speed up processing
setBatchMode(true);

// Create final results table
resultRow = 0;
Table.create("Cytoskeletal Organization - Results");

// ======================================
// MAIN IMAGE ANALYSIS LOOP
// ======================================

for (i = 0; i < validFiles.length; i++) {

    imageName = validFiles[i];
    imagePath = dir + imageName;
    
	baseName = imageName;
	dot = lastIndexOf(baseName, ".");
	if (dot != -1) {
   		baseName = substring(baseName, 0, dot);
	}
    
    print("Analyzing image " + (i+1) + "/" + validFiles.length + ": " + imageName);

    // ======================================
    // CELL MORPHOLOGY ANALYSIS
    // ======================================
    	
   		open(imagePath);
   		run("Set Scale...", "distance=1 known=" + pixelSize + " unit=um global");

        run("8-bit");
         	
        setThreshold(1, 255);
        setOption("BlackBackground", true);
        run("Convert to Mask");
        run("Fill Holes");
        
        run("Set Measurements...", "area mean centroid center perimeter bounding fit shape feret's area_fraction redirect=None decimal=3");
        run("Clear Results");
        run("Analyze Particles...", "clear");

		CellArea = getResult("Area", 0);
		Perimeter = getResult("Perim.", 0);
		MajorAxis = getResult("Major", 0);
		MinorAxis = getResult("Minor", 0);
		AspectRatio = getResult("AR", 0);
		Roundness = getResult("Round", 0);
		Solidity = getResult("Solidity", 0);
		Angle = getResult("Angle", 0);
		
		close("*");
		

    // ======================================
    // FILAMENT THICKNESS MEASUREMENT
    // ======================================
	
	if (measureThickness) {
		
		minBranchThickness_px = minBranchThickness / pixelSize;
		
		open(imagePath);
		
		run("8-bit");
		
		if(IntensityThreshold == 0){
			
			setAutoThreshold("Otsu dark");
			setOption("BlackBackground", true);
			run("Convert to Mask");
			run("Local Thickness (complete process)", "threshold=128");

		} else {
			
			run("Local Thickness (complete process)", "threshold=" + IntensityThreshold);
		
		}
		
		width = getWidth();
		height = getHeight();

		thicknessValues = newArray(width * height);

		sumThickness = 0;
		nPixels = 0;

		for (x = 0; x < width; x++) {
  		  for (y = 0; y < height; y++) {

       	 	value = getPixel(x, y);

        	if (value > minBranchThickness_px) {
        		
        		thicknessValues[nPixels] = value;
            	sumThickness += value;
            	nPixels++;
        		
        		} else {
        			
        		setPixel(x, y, 0);
        			
        		}
    		}
		}

		thicknessValues = Array.trim(thicknessValues, nPixels);
		Array.sort(thicknessValues);

		if (nPixels % 2 == 1) {

    		median_px = thicknessValues[(nPixels-1)/2];

		} else {

    		median_px = (thicknessValues[nPixels/2 - 1] + thicknessValues[nPixels/2]) / 2;
		}

		median_um = median_px * pixelSize;
		
		q1_px = thicknessValues[floor(0.25 * (nPixels - 1))];
		q3_px = thicknessValues[floor(0.75 * (nPixels - 1))];
		q1_um = q1_px * pixelSize;
		q3_um = q3_px * pixelSize;

		meanThickness_px = sumThickness / nPixels;
		meanThickness_um = meanThickness_px * pixelSize;
		
		thicknessMin_px = thicknessMin / pixelSize;
		thicknessMax_px = thicknessMax / pixelSize;
		
		setMinAndMax(thicknessMin_px, thicknessMax_px);
		run("RGB Color");
		
		saveAs("Tiff", outputDir + "Thickness_" + imageName);
		
		close("*");
	}


    // ======================================
    // SKELETON GENERETION
    // ======================================
    
   		open(imagePath);
   		
   		run("8-bit");
   		
   		if (automaticSigma) {
   			
   			sigma = meanThickness_um/2;
			run("Tubeness", "sigma=" + sigma + " use");
			setOption("ScaleConversions", true);
			
		} else {
    
    		run("Tubeness", "sigma=" + sigmaTubeness + " use");
			setOption("ScaleConversions", true);
    
		}
   		
   		run("8-bit");
   		setThreshold(1, 255);
        setOption("BlackBackground", true);
        run("Convert to Mask");
       
       	run("Create Selection");		
		getStatistics(area, mean, min, max, std, histogram);
		CytoskeletonArea = area;
    
    	run("Skeletonize");
    	
    	run("Clear Results");
		roiManager("Reset");
    	
    	run("Analyze Particles...", "size=0-infinity clear add");
    	
    	nROI = roiManager("count");

		for (p = nROI-1; p >= 0; p--) {

  		  	perimeterROI = getResult("Perim.", p);
  		  	AreaROI = getResult("Area", p);

    		if (perimeterROI < 2.1*minBranchLength || AreaROI < 0.05*minBranchLength) {

        	roiManager("Select", p);

        	setForegroundColor(0,0,0);
       		run("Fill");
       		
    		}
		}
    	
		roiManager("Reset");
		run("Remove Overlay");
		run("Select None");

    	saveAs("Tiff", outputDir + "Skeleton_" + imageName);
    	SkeletonPath = outputDir + "Skeleton_" + imageName;
    	
    	if (saveEnhancedSkeleton) {
    		run("Maximum...", "radius=3");
			saveAs("Tiff", outputDir + "Improved_Skeleton_" + imageName);
		}
		
		close("*");
    
    
    // ======================================
    // SKELETON ANALYSIS
    // ======================================    
    	
    	run("Clear Results");
    	
    	open(SkeletonPath);
    	
    	run("Analyze Skeleton (2D/3D)", "prune=none calculate show");
    
    	if (saveIndividualFilaments) {
			
			if (validFiles.length > 1) {

        		singleFilamentsDir = outputDir + "Single_filaments" + File.separator;
        		
        		if (!File.exists(singleFilamentsDir))
            		File.makeDirectory(singleFilamentsDir);
            		
            		selectWindow("Branch information");
					saveAs("Results", singleFilamentsDir + "Single_filament_" + baseName + ".csv");
	
			} else {
				
					selectWindow("Branch information");
					saveAs("Results", outputDir + "Single_filament_" + baseName + ".csv");			
			}
			
		}    
    
    	if (saveGroupsFilaments) {
    		
    		if (validFiles.length > 1) {

        		GroupedFilamentsDir = outputDir + "Grouped_filaments" + File.separator;
        		
        		if (!File.exists(GroupedFilamentsDir))
            		File.makeDirectory(GroupedFilamentsDir);
            		
					selectWindow("Results");
					saveAs("Results", GroupedFilamentsDir + "Grouped_filament_" + baseName + ".csv");	
	
			} else {
				
					selectWindow("Results");
					saveAs("Results", outputDir + "Grouped_filament_" + baseName + ".csv");		
			}
    		
		}
    	
    	AverageBranchLength = 0;
    	totalBranches = 0;
		totalJunctions = 0;
		
		for (a = 0; a < nResults; a++) {

		AverageBranchLength = getResult("Average Branch Length", a);

		if (AverageBranchLength > minBranchLength){
		
    			totalBranches += getResult("# Branches", a);
    			totalJunctions += getResult("# Junctions", a);
			
			}
		}
    	
  		if (saveIndividualFilaments) {
  			TableName = "Single_filament_" + baseName + ".csv";
  		}
  		else {TableName = "Branch information";}
    	
    	size = Table.size(TableName);
    	
    	nBranches = 0;
		sumLength = 0;
		sumLength2 = 0;
		maxLength = 0;
		sumEuclidean = 0;
		nShort = 0;
		nIntermediate = 0;
		nLong = 0;
		sumCos2 = 0;
		sumSin2 = 0;

		for (b = 0; b < size; b++) {

			Length = Table.get("Branch length", b, TableName);
    		Euclidean = Table.get("Euclidean distance", b, TableName);
    		Intensity = Table.get("average intensity", b, TableName);
    		V1x = Table.get("V1 x", b, TableName);
			V1y = Table.get("V1 y", b, TableName);
			V2x = Table.get("V2 x", b, TableName);
			V2y = Table.get("V2 y", b, TableName);

   			if (Length > minBranchLength &&
        		Euclidean != 0 &&
       			Intensity != 0) {

   					nBranches++;

    				sumLength += Length;
    				
    				sumLength2 += Length * Length;
    				
    				sumEuclidean += Euclidean;
					
					if (Length > maxLength){
    						maxLength = Length;
    				}
       			 
       				if (Length < ShortIntermediate){nShort++;}
					else if (Length < IntermediateLong){nIntermediate++;}
					else{nLong++;}
       			 
       			  	dx = V2x - V1x;
    			  	dy = V2y - V1y;
					AngleWithXaxis = atan2(dy, dx);
    				cos2 = cos(2 * AngleWithXaxis);
    				sin2 = sin(2 * AngleWithXaxis);
    				sumCos2 += cos2;
    				sumSin2 += sin2;
       			 
       			}
		}
   		
   		meanLength = sumLength/nBranches;
    	
    	variance = (sumLength2 - sumLength*sumLength/nBranches)/(nBranches-1);
		stdLength = sqrt(variance);
    	
    	meanTortuosity = sumLength / sumEuclidean;
    	
    	meanCos2 = sumCos2 / nBranches;
		meanSin2 = sumSin2 / nBranches;
		orientationOrder = sqrt(meanCos2*meanCos2 + meanSin2*meanSin2);
		
		meanOrientation = 0.5 * atan2(meanSin2, meanCos2);
		meanOrientationDegrees = meanOrientation*180/PI;
    	
    	selectWindow(TableName);
    	run("Close");
    	selectWindow("Results");
    	run("Close");
    	
    
    // ======================================
    // BRANCHING ANALYSIS
    // ======================================  
    
 		selectImage("Longest shortest paths");
    	
    	w = getWidth();
		h = getHeight();

		for (y = 1; y < h-1; y++) {
    		for (x = 1; x < w-1; x++) {
    
			if (getPixel(x,y) == 255) {
    			connected = false;

    		for (dy=-1; dy<=1; dy++) {
        		for (dx=-1; dx<=1; dx++) {

            if (dx==0 && dy==0)
                continue;

            if (getPixel(x+dx, y+dy) == 255)
                connected = true;
        	}
    	}

    		if (!connected)
        		setPixel(x, y, 96);
			}
   		}
	}
    	    	
        if (saveBranchMap || saveEnhancedBranchMap) {
			
			run("Duplicate...", "title=temp");
			run("RGB Color");
			
			for (y=0; y<h; y++) {
  			  for (x=0; x<w; x++) {

        		value = getPixel(x,y);

        		r = (value & 0xff0000) >> 16;
      			g = (value & 0x00ff00) >> 8;
       			b = value & 0x0000ff;

	        if (r==255 && g==255 && b==255) {
    	        setPixel(x,y,0xffff00);
        		}

	        if (r==195 && g==0 && b==93) {
    	        setPixel(x,y,0x00ffff);
        		}
    		}
		}
		
		if (saveBranchMap) {
			saveAs("Tiff", outputDir + "Branching_" + imageName);
			}	
		
		if (saveEnhancedBranchMap) {
			
			run("Maximum...", "radius=3");	
			
			for (y=0; y<h; y++) {
  			  for (x=0; x<w; x++) {

        		value2 = getPixel(x,y);

        		r = (value2 & 0xff0000) >> 16;
      			g = (value2 & 0x00ff00) >> 8;
       			b = value2 & 0x0000ff;

	        if (r==255 && g==255 && b==255) {
    	        setPixel(x,y,0x00ffff);
        		}
    		}
		}
				
			saveAs("Tiff", outputDir + "Improved_Branching_" + imageName);
			}	
		
		close();	
		
		}	
		
		MainBranches = 0;
		SideBranches = 0;

		for (y = 1; y < h-1; y++) {
    		for (x = 1; x < w-1; x++) {

				if (getPixel(x,y) == 96) {
					MainBranches++;
				}


				if (getPixel(x,y) == 255) {
					SideBranches++;
				}
    		}
		}	

    	close("*");
    	
    // ======================================
    // ORIENTATION ANALYSIS
    // ====================================== 
    	
    	open(SkeletonPath);
    	
    	run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
		count_total = Table.get("Count", 0, "Summary");

		gradient = -1;

		for (g = 0; g < methods.length; g++) {
   	 		if (gradientMethod == methods[g]) {
        		gradient = g;
        		break;
    		}
		}

		run("OrientationJ Analysis", "tensor=" + orientationSigma + " gradient=" + gradient + " color-survey=on hsb=on hue=Orientation sat=Constant bri=Original-Image radian=on ");

		if (saveOrientationMap) {
			saveAs("Tiff", outputDir + "Orientation_" + imageName);
		}	
		
		if (saveEnhancedOrientationMap) {
			run("Duplicate...", "title=Temp");
			run("Maximum...", "radius=3");	
			saveAs("Tiff", outputDir + "Improved_Orientation_" + imageName);
			close("Temp");
		}	

		run("HSB Stack");
		run("Stack to Images");

		count_segment = 0;
		counts = newArray(256);

		for (b = 0; b < 256; b++) {

			AngleInterval = floor(directionChangeThreshold*256/180);

    		min = b;
    		max = (b + AngleInterval - 1) % 256;
    		
			selectImage("Hue");
			run("Duplicate...", "title=Hue_temp");
			selectImage("Hue_temp");

    	if (min <= max) {

			if (b == 0){
			        
			setThreshold(0.1, max);
        	run("Convert to Mask");
        
        	run("Clear Results");
			run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
			count_segment = Table.get("Count", b + 1, "Summary");
		}

		else {

        	setThreshold(min, max);
        	run("Convert to Mask");
        
        	run("Clear Results");
			run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
			count_segment = Table.get("Count", b + 1, "Summary");
		}

    } else {

        	setThreshold(min, 255);
        	run("Convert to Mask");
        	rename("part1");

			selectImage("Hue");
			run("Duplicate...", "title=part2");
			selectImage("part2");

        	setThreshold(0.1, max);
        	run("Convert to Mask");

        	imageCalculator("OR create", "part1", "part2");
        
        	run("Clear Results");
			run("Analyze Particles...", "size=3-Infinity pixel clear summarize");
			count_segment = Table.get("Count", b + 1, "Summary");

        	close("part1");
        	close("part2");

    	}

		counts[b] = count_segment;
	
    	close();
    	
	}

	selectWindow("Summary");
    run("Close");

	sumOrientation = newArray(AngleInterval);
	nBins = floor(256 / AngleInterval);

	for (shift = 0; shift < AngleInterval; shift++) {

    	totalOrientation = 0;

    	for (q = 0; q < nBins; q++) {

        	index = shift + q * AngleInterval;

        	totalOrientation += counts[index];

    	}

    	sumOrientation[shift] = totalOrientation;
		
	}

	meanNofSegments = 0;

	for (m = 0; m < AngleInterval; m++) {
    	meanNofSegments += sumOrientation[m];
	}

	meanNofSegments = meanNofSegments / AngleInterval;
    NofSegments = meanNofSegments - count_total;
	
	close("*");

    // ======================================
    // RESULTS
    // ======================================
    
    CellCytoskeletonArea = CytoskeletonArea / CellArea * 100;
    
    JunctionOnBranches = totalJunctions/nBranches;
    
    FilamentDensity = nBranches/CellArea * 100;
    
    ShortFilamentDensity = nShort/CellArea * 100;
    IntermediateFilamentDensity = nIntermediate/CellArea * 100;
    LongFilamentDensity = nLong/CellArea * 100;
    
    FractionBranchedFilaments = (SideBranches/(SideBranches+MainBranches))*100;
    
    N_change_in_direction = NofSegments/sumLength;	
    
    
	AreaHeader = "Cell Area (" + micro + "m" + sup2 + ")";
    MeanThicknessHeader = "Mean Thickness (" + micro + "m" + ")";
    Q1ThicknessHeader = "Q1 Thickness (" + micro + "m" + ")";
    MedianThicknessHeader = "Median Thickness (" + micro + "m" + ")";
    Q3ThicknessHeader = "Q3 Thickness (" + micro + "m" + ")";
    CytoskeletonAreaHeader = "Cytoskeleton Area (" + micro + "m" + sup2 + ")";
    AcskAcellHeader = "Acsk/Acell (%)";
    LengthHeader = "Mean filament length (" + micro + "m)";
	StdLengthHeader = "Std Dev filament length (" + micro + "m)";
	TotalLengthHeader = "Total filament length (" + micro + "m)";
	MaxLengthHeader = "Maximum filament length (" + micro + "m)";
	FilamentDensityHeader = "Filament Density per 100 " + micro + "m" + sup2;
    ShortHeader = "N° of Short Filaments (<" + ShortIntermediate + " " + micro + "m)";
    IntermediateHeader = "N° of Intermediate Filaments (" + ShortIntermediate + " - " + IntermediateLong + " " + micro + "m)"; 
    LongHeader = "N° of Long Filaments (>" + IntermediateLong + " " + micro + "m)";
	ShortDensityHeader = "Short Filament per 100 " + micro + "m" + sup2;
	IntermediateDensityHeader = "Intermediate Filament per 100 " + micro + "m" + sup2;
	LongDensityHeader = "Long Filaments per 100 " + micro + "m" + sup2;
	DirectionChangeHeader = "Direction Changes per " + micro + "m"; 
    
    Table.set("Image", resultRow, baseName, "Cytoskeletal Organization - Results");

    Table.set(AreaHeader, resultRow, CellArea, "Cytoskeletal Organization - Results");
	Table.set("Cell Roundness", resultRow, Roundness, "Cytoskeletal Organization - Results");
	Table.set("Cell Solidity", resultRow, Solidity, "Cytoskeletal Organization - Results");	
    
	if(measureThickness){
		Table.set(MeanThicknessHeader, resultRow, meanThickness_um, "Cytoskeletal Organization - Results");
		Table.set(Q1ThicknessHeader, resultRow, q1_um, "Cytoskeletal Organization - Results");
		Table.set(MedianThicknessHeader, resultRow, median_um, "Cytoskeletal Organization - Results");
		Table.set(Q3ThicknessHeader, resultRow, q3_um, "Cytoskeletal Organization - Results");
	}

	Table.set(CytoskeletonAreaHeader, resultRow, CytoskeletonArea, "Cytoskeletal Organization - Results");
	Table.set(AcskAcellHeader, resultRow, CellCytoskeletonArea, "Cytoskeletal Organization - Results");
	
	Table.set("N° of Filaments", resultRow, nBranches, "Cytoskeletal Organization - Results");
	Table.set("N° of Junctions", resultRow, totalJunctions, "Cytoskeletal Organization - Results");
	Table.set("Junction/Filaments", resultRow, JunctionOnBranches, "Cytoskeletal Organization - Results");
	Table.set(LengthHeader, resultRow, meanLength, "Cytoskeletal Organization - Results");
	Table.set(StdLengthHeader, resultRow, stdLength, "Cytoskeletal Organization - Results");
	Table.set(TotalLengthHeader, resultRow, sumLength, "Cytoskeletal Organization - Results");
	Table.set(MaxLengthHeader, resultRow, maxLength, "Cytoskeletal Organization - Results");
	Table.set("Mean Tortuosity", resultRow, meanTortuosity, "Cytoskeletal Organization - Results");
	Table.set(FilamentDensityHeader, resultRow, FilamentDensity, "Cytoskeletal Organization - Results");
	Table.set(ShortHeader, resultRow, nShort, "Cytoskeletal Organization - Results");
	Table.set(IntermediateHeader, resultRow, nIntermediate, "Cytoskeletal Organization - Results");
	Table.set(LongHeader, resultRow, nLong, "Cytoskeletal Organization - Results");
	Table.set(ShortDensityHeader, resultRow, ShortFilamentDensity, "Cytoskeletal Organization - Results");
	Table.set(IntermediateDensityHeader, resultRow, IntermediateFilamentDensity, "Cytoskeletal Organization - Results");
	Table.set(LongDensityHeader, resultRow, LongFilamentDensity, "Cytoskeletal Organization - Results");
	Table.set("Mean Orientation (Degrees)", resultRow, meanOrientationDegrees, "Cytoskeletal Organization - Results");
	Table.set("Orientation Order Parameter", resultRow, orientationOrder, "Cytoskeletal Organization - Results");
	Table.set("Fraction of Branched Filaments (%)", resultRow, FractionBranchedFilaments, "Cytoskeletal Organization - Results");
	Table.set(DirectionChangeHeader, resultRow, N_change_in_direction, "Cytoskeletal Organization - Results");
	
	Table.update("Cytoskeletal Organization - Results");
	
	resultRow++;

}

selectWindow("Cytoskeletal Organization - Results");
saveAs("Results", outputDir + "Results.csv");

setBatchMode(false);
print("Analysis Completed");