function CellType = PreProcess_Step4_CellType_Lung_v2(channels)

ttf1  = channels.ttf1 ;
panck = channels.panck ;

cd45  = channels.cd45 ;

% cd19  = channels.cd19 ;
cd20  = channels.cd20 ;

cd3   = channels.cd3 ;
cd4   = channels.cd4  ;
cd8   = channels.cd8 ;
foxp  = channels.foxp ;

cd14  = channels.cd14 ;
cd16  = channels.cd16 ;
cd163 = channels.cd163 ;
cd68  = channels.cd68 ;

asma  = channels.asma ;
vim   = channels.vim ;

% conflict resolution threshold
CellType.thresh = 100;
CellType.index = [];
CellType.codes  = [];
CellType.layer = [];
CellType.names  = {};
CellType.layerjump = 10;  % either 10 or 100 usually

CellType.NameList = {{'Immune','Epithelial','Stroma','Other'} ...
                    ,{'Lymphoid','Myeloid'} ...
                    ,{'T','B','TAM'} ...
                    ,{'T reg','T helper','T cytotox','TAM CD68+','TAM CD163+'} ...
                    };

CellType.Classes =  { ...
              {'Immune',           1, 1     ,{cd45,cd20,cd3,cd4,cd8,cd14,cd16,cd68,cd163}} ...   
             ,{'Epithelial',       1, 2     ,{ttf1,panck}}  ...
             ,{'Stroma',           1, 3     ,{[asma -cd45 -cd20 -cd3 -cd4 -cd8 -cd14 -cd16 -cd68 -cd163], [vim -cd45 -cd20 -cd3 -cd4 -cd8 -cd14 -cd16 -cd68 -cd163]}}  ...
             ,{'Other',            1, 4     ,{[-ttf1 -panck -asma -vim -cd45 -cd20 -cd3 -cd4 -cd8 -cd14 -cd16 -cd68 -cd163]}} ... 
             ,{'Lymphoid',         2, 11    ,{cd20,cd3,[cd3 cd4],cd8}}...
             ,{'Myeloid',          2, 12    ,{cd14,cd16,cd68,cd163}}        ...
             ,{'T',                3, 111   ,{cd3,[cd3 cd4],cd8}}     ... 
             ,{'B',                3, 112   ,{cd20}}               ... 
             ,{'TAM',              3, 121   ,{cd14,cd68,cd163}}  ... 
             ,{'T reg',            4, 1111  ,{[cd4 foxp],[cd8 foxp]}}             ...
             ,{'T helper',         4, 1112  ,{[cd4 -foxp]}}         ... 
             ,{'T cytotox',        4, 1113  ,{[cd8 -foxp]}}           ... 
             ,{'TAM CD68+',        4, 1211  ,{cd68}}           ... 
             ,{'TAM CD163+',       4, 1212  ,{cd163}}           ... 
             };

% ABlist is an ordered list of all the channels (antibodies) used for cell typing
CellType.ABlist = cellfun(@(x) x{4},CellType.Classes,'un',0);
CellType.ABlist = unique(cell2mat(cellfun(@(x) abs([x{:}]),CellType.ABlist,'un',0)));
% use this list to specify which markers are Nuclear vs Cytoplasmic
    % 1 for nuclear (default), 0 for cytoplasmic
CellType.NvC = zeros(size(CellType.ABlist)); 
    % for manual override of all markers
CellType.NvC(CellType.ABlist==ttf1) = 1;
CellType.NvC(CellType.ABlist==foxp) = 1;
CellType.NvC(CellType.ABlist==asma) = 1;

for i = 1:length(CellType.Classes)
    CellType.index = [CellType.index; i];
    CellType.codes  = [CellType.codes;  CellType.Classes{i}{3}];
    CellType.layer = [CellType.layer; CellType.Classes{i}{2}];
    CellType.names  = [CellType.names;  CellType.Classes{i}{1}];
end