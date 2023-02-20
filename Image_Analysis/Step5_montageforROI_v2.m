function Step5_montageforROI_v2(filename,options)
% This program saves the smallest montage made by ASHLAR for ROI drawing on
% ImageJ

for i = 1:length(filename.folders.fols)
    tic
    fold = filename.folders.fols{i};
    disp(fold) 
    
    %Making folders
    outputfolder = [filename.folders.main filename.folders.output filename.folders.montagelowres];
    addpath(filename.folders.main)
    mkdir(outputfolder)
  
    filenamesv = [outputfolder fold '_montage.tif'];
    omifile = [ filename.folders.main filename.folders.ometiff filename.folders.fols{i} filename.ashlarsuffix];
    try
        ominfo = imfinfo(omifile); 
    catch
        disp(['Ometiff file ' filename.folders.fols{i} '.ome.tif not found'])
        continue
    end
    slices = length(ominfo); %Total number of slices in the omi.tif file 
    
    smallslice = filename.maxround*(options.pyramidlevel-1)+1; 
    
    imwrite(DAPImontage, filenamesv) 
    toc 
end 

