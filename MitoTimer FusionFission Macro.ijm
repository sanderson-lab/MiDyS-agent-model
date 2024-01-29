green_dir = getDirectory("Green Directory");
red_dir = getDirectory("Red Directory");
source_dir = getDirectory("Source Directory");

results_dir = source_dir + File.separator + "Results";
File.makeDirectory(results_dir);

green_list = getFileList(green_dir);
red_list = getFileList(red_dir);

green_seg = getDirectory("Green Seg Directory");
red_seg = getDirectory("Red Seg Directory");

setBatchMode(true);

green_seg_list = getFileList(green_seg);
red_seg_list = getFileList(red_seg);

run("Set Measurements...", "integrated display redirect=None decimal=3");

mito_dir = source_dir + File.separator + "Mitochondria";
File.makeDirectory(mito_dir);

for (j=0; j<green_seg_list.length; j++)
{
	open(green_seg + File.separator + green_seg_list[j]);
	img = getTitleStripExtension();
	run("8-bit");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Despeckle");
	rename("green" + j);
	run("Measure");
	open(red_seg + File.separator + red_seg_list[j]);
	run("8-bit");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Despeckle");
	rename("red" + j);
	run("Measure");
	imageCalculator("Add create", "green" + j, "red" + j);
	run("Canvas Size...", "width=1360 height=1030 position=Top-Left");
	run("Morphological Filters", "operation=Erosion element=Disk radius=0.5");
	saveAs("Tiff", mito_dir + File.separator + img + "Mitochondria");
	close("*");
}

saveAs("Results", results_dir + File.separator + "IntDenResults.csv");
close("Results");

mito_list = getFileList(mito_dir);

fission_dir = source_dir + File.separator + "Fission Events";
File.makeDirectory(fission_dir);

fission_roi = source_dir + File.separator + "Fission ROIs";
File.makeDirectory(fission_roi);

for (k=0; k< (mito_list.length-1); k++)
{
	print("Processing Slice: " + k);
    open(mito_dir + File.separator + mito_list[k]);
    img = getTitleStripExtension();
    rename("topslice");
	run("Analyze Particles...", "size=0.005-Infinity show=Nothing add");
	roiManager("Set Color", "red");
	roiManager("Set Line Width", 1);
	open(mito_dir + File.separator + mito_list[k+1]);
	rename("bottomslice");
	roiManager("Show None");
	nROI = roiManager("count");
	fission_events = nROI;
	non_events = newArray();
	for (m=0; m < nROI; m++)
	{
		roiManager("select", m);
		run("Analyze Particles...", "size=0.001-Infinity summarize");
		IJ.renameResults("Summary","Results");
		fis = getResult("Count", 0);
		if (fis < 2) {
			non_events = Array.concat(non_events, m);
			fission_events--;
		}
		close("Results");
		roiManager("deselect");
	}
	roiManager("select", non_events);
	roiManager("delete");
	roiManager("Show All without labels");
	roiManager("Save", fission_roi + File.separator + img + "fission_ROIs.zip");
	run("Flatten");
	saveAs("Tiff", fission_dir + File.separator + img + "Slice");
	close("*");
	roiManager("delete");
	print("Fission Events: " + fission_events);
	fission_rate = fission_events / nROI * 100;
	print("Fission Rate: " + fission_rate + "%");
}

fis_roi_list = getFileList(fission_roi);

for (l=0; l<fis_roi_list.length; l++)
{
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=3");
	open(mito_dir + File.separator + mito_list[l]);
	roiManager("Open", fission_roi + File.separator + fis_roi_list[l]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "FissionShapes.csv");
close("Results");

for (l=0; l<fis_roi_list.length; l++)
{
	run("Set Measurements...", "mean integrated median display redirect=None decimal=3");
	open(green_dir + File.separator + green_list[l]);
	roiManager("Open", fission_roi + File.separator + fis_roi_list[l]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "FissionGreen.csv");
close("Results");

for (l=0; l<fis_roi_list.length; l++)
{
	open(red_dir + File.separator + red_list[l]);
	roiManager("Open", fission_roi + File.separator + fis_roi_list[l]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "FissionRed.csv");
close("Results");

total_roi = source_dir + File.separator + "Total ROIs";
File.makeDirectory(total_roi);

for (o=0; o<mito_list.length; o++)
{
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=3");
	open(mito_dir + File.separator + mito_list[o]);
	img = getTitleStripExtension();
	run("Analyze Particles...", "size=0.005-Infinity display add");
	roiManager("Save", total_roi + File.separator + img + "total_ROIs.zip");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "TotalShapes.csv");
close("Results");

total_roi_list = getFileList(total_roi);

for (o=0; o<green_list.length; o++)
{
	run("Set Measurements...", "mean integrated median display redirect=None decimal=3");
	open(green_dir + File.separator + green_list[o]);
	getTitleStripExtension();
	roiManager("Open", total_roi + File.separator + total_roi_list[o]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "TotalGreen.csv");
close("Results");

for (o=0; o<red_list.length; o++)
{
	run("Set Measurements...", "mean integrated median display redirect=None decimal=3");
	open(red_dir + File.separator + red_list[o]);
	getTitleStripExtension();
	roiManager("Open", total_roi + File.separator + total_roi_list[o]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "TotalRed.csv");
close("Results");

fusion_dir = source_dir + File.separator + "Fusion Events";
File.makeDirectory(fusion_dir);

fusion_roi = source_dir + File.separator + "Fusion ROIs";
File.makeDirectory(fusion_roi);

fusion_mask = source_dir + File.separator + "Fusion Masks";
File.makeDirectory(fusion_mask);

fusion_original_roi = source_dir + File.separator + "Fusion Original ROIs";
File.makeDirectory(fusion_original_roi);

for (k=0; k< (mito_list.length-1); k++)
{
	print("Processing Slice: " + k);
    open(mito_dir + File.separator + mito_list[k+1]);
    img = getTitleStripExtension();
    rename("bottomslice");
	run("Analyze Particles...", "size=0.005-Infinity show=Nothing add");
	roiManager("Set Color", "red");
	roiManager("Set Line Width", 1);
	close("*");
	open(mito_dir + File.separator + mito_list[k]);
	rename("topslice");
	roiManager("Show None");
	nROI = roiManager("count");
	fusion_events = nROI;
	non_events2 = newArray();
	for (m=0; m < nROI; m++)
	{
		roiManager("select", m);
		run("Analyze Particles...", "size=0.001-Infinity summarize");
		IJ.renameResults("Summary","Results");
		fus = getResult("Count", 0);
		if (fus < 2) {
			non_events2 = Array.concat(non_events2, m);
			fusion_events--;
		}
		close("Results");
		roiManager("deselect");
	}
	roiManager("select", non_events2);
	roiManager("delete");
	roiManager("Show All without labels");
	roiManager("Save", fusion_roi + File.separator + img + "fusion_ROIs.zip");
	run("Flatten");
	saveAs("Tiff", fusion_dir + File.separator + img + "Slice");
	close("*");
	print("Fusion Events: " + fusion_events);
	fusion_rate = fusion_events / nROI * 100;
	print("Fusion Rate: " + fusion_rate + "%");
	open(mito_dir + File.separator + mito_list[k]);
	rename("mito");
	nROI = roiManager("count");
	for (m=0; m < nROI; m++)
	{
		selectWindow("mito");
		roiManager("select", m);
		run("Analyze Particles...", "  show=Masks");
	}
	close("mito");
	roiManager("deselect");
	roiManager("delete");
	run("Images to Stack", "name=Stack title=[] use");
	run("Z Project...", "projection=[Sum Slices]");
	run("8-bit");
	saveAs("Tiff", fusion_mask + File.separator + img + "Slice");
	close("Stack");
	open(mito_dir + File.separator + mito_list[k]);
	run("Analyze Particles...", "size=0.00-Infinity show=Nothing add");
	close();
	selectWindow(img + "Slice.tif");
	nROI = roiManager("count");
	run("Set Measurements...", "mean display redirect=None decimal=3");
	roiManager("Measure");
	non_events3 = newArray();
	for (m=0; m < nROI; m++)
	{
		particle = getResult("Mean", m);
		if (particle == 0) {
			non_events3 = Array.concat(non_events3, m);
		}
	}
	close("Results");
	roiManager("select", non_events3);
	roiManager("delete");
	roiManager("Save", fusion_original_roi + File.separator + img + "fusion_ROIs.zip");
}

fus_original_roi_list = getFileList(fusion_original_roi);

for (l=0; l<fus_original_roi_list.length; l++)
{
	run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=3");
	open(mito_dir + File.separator + mito_list[l]);
	roiManager("Open", fusion_original_roi + File.separator + fus_original_roi_list[l]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "FusionShapes.csv");
close("Results");

for (l=0; l<fus_original_roi_list.length; l++)
{
	run("Set Measurements...", "mean integrated median display redirect=None decimal=3");
	open(green_dir + File.separator + green_list[l]);
	roiManager("Open", fusion_original_roi + File.separator + fus_original_roi_list[l]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "FusionGreen.csv");
close("Results");

for (l=0; l<fus_original_roi_list.length; l++)
{
	open(red_dir + File.separator + red_list[l]);
	roiManager("Open", fusion_original_roi + File.separator + fus_original_roi_list[l]);
	roiManager("measure");
	roiManager("delete");
	close("*");
}
saveAs("Results", results_dir + File.separator + "FusionRed.csv");
close("Results");


function getTitleStripExtension() {
  t = getTitle();
  t = replace(t, ".tif", "-");
  return t;
}