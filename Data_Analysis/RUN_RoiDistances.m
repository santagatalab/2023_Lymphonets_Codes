clear all
MasterFolder = ['Y:\sorger\data\RareCyte\Giorgio\Mouse KP Lung CycIF Master Folder\'];

basefolder = [MasterFolder '2019-09 KP LucOS Cxcl10\'];
analfolder = [basefolder 'ANALYSIS\DATA ANALYSIS MAY2021\'];
resufolder = 'Results_singletissue\'; 
img_analfolder = [basefolder 'ANALYSIS\IMAGE ANALYSIS\'];
date = '20210526';

codedir = [MasterFolder 'common functions'];
addpath(codedir)

load([analfolder resufolder 'Results_Morp_' date '.mat'])
load([analfolder resufolder 'Results_Settings_' date '.mat'])

filename.basefolder = basefolder;
filename.analfolder = analfolder;
filename.resufolder = resufolder;
options.date = date;
options.figOpt = 0;
% save([analfolder resufolder 'Results_Settings_' date '.mat'],'filename','options')

%% calculate the distance of cells from the border of the closest tumor and compare it to the distance from the closest blood vessel

% strategy:
% - load mask from tumors and blood vessels
% 1) measure the distance from blood vessel
% 2) measure the distance from a tumor
% 3) measure the distance from the border of a tumor

DistResults.Tumor = zeros(size(MorpResults.Indexes,1),3); 
% col (1) tumor #
% col (2) distance from tumor - 0 if inside
% col (3) distance from tumor boundary
DistResults.Blood = zeros(size(MorpResults.Indexes,1),1); 
% col (1) distance from blood vessel - 0 if inside
% col (2) distance from tumor boundary

for tis = 1:max(MorpResults.Indexes)
    tis
    % load masks for tumor and blood
    basename = [img_analfolder filename.roifolder filename.folders{tis}(5:end) '_'];
    TumorMask = imread([basename options.tifnames.Path{1}]);
    BloodMask = imread([basename options.tifnames.Path{2}]);
    TumorBoundMask = bwperim(TumorMask);
    BloodBoundMask = bwperim(BloodMask);
    
    % calcuate the distance
    [TumorDist, idx] = bwdist(TumorMask);
    BloodDist = bwdist(BloodMask);
    TumorBoundDist = bwdist(TumorBoundMask);
    BloodBoundDist = bwdist(BloodBoundMask);
    
    % get cell index
    vect_cell_ind = MorpResults.Indexes==tis;
    X_scale = max([MorpResults.X(vect_cell_ind)/(2^(options.ROI.pyramidlevel+1)) ones(sum(vect_cell_ind),1)],[],2);
    Y_scale = max([MorpResults.Y(vect_cell_ind)/(2^(options.ROI.pyramidlevel+1)) ones(sum(vect_cell_ind),1)],[],2);
    cell_index = sub2ind(size(TumorMask),Y_scale,X_scale);
    
    % calculate distance matrices
    DistResults.Tumor(vect_cell_ind,1) = TumorMask(idx(cell_index));
    DistResults.Tumor(vect_cell_ind,2) = TumorDist(cell_index);
    DistResults.Tumor(vect_cell_ind,3) = TumorBoundDist(cell_index);
    
    DistResults.Blood(vect_cell_ind,1) = BloodDist(cell_index);
    DistResults.Blood(vect_cell_ind,2) = BloodBoundDist(cell_index);
    
end
% make the distance from boundary negative if the cell is inside the tumor
DistResults.Tumor(DistResults.Tumor(:,2)==0,3) = -DistResults.Tumor(DistResults.Tumor(:,2)==0,3);
DistResults.Blood(DistResults.Blood(:,1)==0,2) = -DistResults.Blood(DistResults.Blood(:,1)==0,2);


disp('Distance calculation done')
%%% SAVE
save([filename.analfolder filename.resufolder 'Results_RoiDist_' options.date '.mat'],'DistResults')
