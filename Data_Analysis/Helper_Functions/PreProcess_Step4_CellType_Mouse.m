function CellType = PreProcess_Step4_CellType(channels)

ttf1 = channels.ttf1 ;
cd45r = channels.cd45r ;
cd45 = channels.cd45 ;
foxp = channels.foxp ;
cd4  = channels.cd4  ;
cd8  = channels.cd8 ;
cd103 = channels.cd103 ;
cd11c = channels.cd11c ;
cd11b = channels.cd11b ;
nkp46 = channels.nkp46 ;
gzmb = channels.gzmb ;
perf = channels.perf ;
ly6g = channels.ly6g ;

% conflict resolution threshold
CellType.thresh = 100;
CellType.index = [];
CellType.codes  = [];
CellType.layer = [];
CellType.names  = {};
CellType.layerjump = 10;  % either 10 or 100 usually

CellType.NameList = {{'Immune','Epithelial','Other'} ...
                    ,{'Lymphoid','Myeloid'} ...
                    ,{'T','B','NK_L','Alveolar MAC','DC','NK_M','Neutrophil','TAM'} ...
                    ,{'T reg','T helper','T cytotox'} ...
                    };

CellType.Classes =  { ...
              {'Immune',           1, 1     ,{cd45,cd45r,cd4,cd8,cd11b,cd11c}} ...   
             ,{'Epithelial',       1, 2     ,{ttf1}}  ...
             ,{'Other',            1, 3     ,{[-ttf1 -cd45 -cd45r -cd4 -cd8 -cd11b -cd11c -cd103]}} ... 
             ,{'Lymphoid',         2, 11    ,{cd4,cd8,[nkp46 -cd103],cd45r}}...
             ,{'Myeloid',          2, 12    ,{cd11b,cd11c}}        ...
             ,{'T',                3, 111   ,{cd4,cd8}}     ... 
             ,{'B',                3, 112   ,{cd45r}}               ... 
             ,{'NK_L',             3, 113   ,{[nkp46 -cd103]}} ...
             ,{'Alveolar MAC',     3, 121   ,{[cd11c -cd103 -cd4 -cd8]}}  ... 
             ,{'DC',               3, 122   ,{[cd11c cd103],[cd11c cd8]}} ...    
             ,{'NK_M',             3, 123   ,{[cd11b gzmb],[cd11b perf],[nkp46 -cd103]}} ...
             ,{'Neutrophil',       3, 124   ,{[cd11b ly6g]}}              ... 
             ,{'TAM',              3, 125   ,{[cd11b -cd11c -ly6g -gzmb -perf]}}  ... 
             ,{'T reg',            4, 1111  ,{[cd4 foxp],[cd8 foxp]}}             ...
             ,{'T helper',         4, 1112  ,{[cd4 -foxp]}}         ... 
             ,{'T cytotox',        4, 1113  ,{[cd8 -foxp]}}           ... 
              };

% ABlist is an ordered list of all the channels (antibodies) used for cell typing
CellType.ABlist = cellfun(@(x) x{4},CellType.Classes,'un',0);
CellType.ABlist = unique(cell2mat(cellfun(@(x) abs([x{:}]),CellType.ABlist,'un',0)));
% use this list to specify which markers are Nuclear vs Cytoplasmic
    % 1 for nuclear (default), 0 for cytoplasmic
CellType.NvC = zeros(size(CellType.ABlist)); 
    % for manual override of all markers
if ttf1 > 0
    CellType.NvC(CellType.ABlist==ttf1) = 1;
end
if gzmb > 0
    CellType.NvC(CellType.ABlist==gzmb) = 1;
end
if perf > 0
    CellType.NvC(CellType.ABlist==perf) = 1;
end
if foxp > 0
    CellType.NvC(CellType.ABlist==foxp) = 1;
end
if ly6g > 0
    CellType.NvC(CellType.ABlist==ly6g) = 1;
end

for i = 1:length(CellType.Classes)
    CellType.index = [CellType.index; i];
    CellType.codes  = [CellType.codes;  CellType.Classes{i}{3}];
    CellType.layer = [CellType.layer; CellType.Classes{i}{2}];
    CellType.names  = [CellType.names;  CellType.Classes{i}{1}];
end