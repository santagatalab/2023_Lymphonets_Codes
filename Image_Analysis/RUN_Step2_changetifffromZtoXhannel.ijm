DIR ="Z:/sorger/data/IN_Cell_Analyzer_6000/Claire/Mouse_Lung_Megan_CXCL10_Experiment_09-2019/ANALYSIS/Ilastik Segmentation/ANALYSIS/"; 

myType =newArray("CroppedData","FullStacks"); //"_forIlastik",
suffix = ".tif";

for (j=0; j< myType.length; j++) {
		myDIRin = DIR + myType[j];  //+ myDIR[k] 
		print(myDIRin);
		list = getFileList(myDIRin);
		for (i=0; i < list.length; i++) {
			if (endsWith(list[i], suffix)) {
		// Make sure to change channel below 
			open(myDIRin+"\\"+list[i]);
			run("Properties...", "channels=8 slices=1 frames=1 unit=pixel pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");
			run("Save");	
			run("Close");
			}
		}
}

//for (k=0; k< myDIR.length; k++) { }
//myDIR = newArray("MLB_267_6','MLB_268_2","MLB_269_5_ROI1","MLB_269_5_ROI2","MLB_270_1","MLB_271_3","MLB_272_2_ROI1","MLB_272_2_ROI2","MLB_273_4_ROI1","MLB_273_4_ROI2","MLB_274_3","MLB_275_1","MLB_276_2","MLB_277_1","MLB_278_4_ROI1");
//                 //'MLB_278_4_ROI2', 'MLB_279_4_ROI1', 'MLB_279_4_ROI2','MLB_280_2_ROI1', 'MLB_280_2_ROI2', 'MLB_281_3_ROI1',  'MLB_281_3_ROI2', 'MLB_282_4_ROI1', 'MLB_282_4_ROI2','MLB_283_1',      'MLB_284_2',      'MLB_285_6',  'MLB_266_1');

