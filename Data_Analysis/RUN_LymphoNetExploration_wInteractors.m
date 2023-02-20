clear all
MasterFolder = ['Y:\sorger\data\RareCyte\Giorgio\Mouse KP Lung CycIF Master Folder\'];

basefolder = [MasterFolder '2019-09 KP LucOS Cxcl10\'];
analfolder = [basefolder 'ANALYSIS\DATA ANALYSIS MAY2021\'];
resufolder = 'Results_singletissue\'; 
date = '20210526';

codedir = [MasterFolder 'common functions'];
addpath(codedir)

load([analfolder resufolder 'Results_Morp_' date '.mat'])
load([analfolder resufolder 'Results_ROI_20200408.mat'])
load([analfolder resufolder 'Results_Norm_' date '.mat']);
load([analfolder resufolder 'Results_Filt_' date '.mat']);
load([analfolder resufolder 'Results_CellType_' date '.mat'])
load([analfolder resufolder 'Results_Settings_' date '.mat'])

filename.basefolder = basefolder;
filename.analfolder = analfolder;
filename.resufolder = resufolder;
options.date = date;
options.figOpt = 0;
save([analfolder resufolder 'Results_Settings_' date '.mat'],'filename','options')


%% Settings

distance_cutoff = [25 50];

for d = 1:length(distance_cutoff)
    LymphoNets.thresh.dist = distance_cutoff(d) ; 
    LymphoNets.thresh.netsize = 5;
    % Cell Types to look for    ***will also be used for the network filter***
    LymphoNets.thresh.CellTypes = ... Row 1: CellType to look for, Row 2: CellType label
        {'B','T','Treg','Thelp','Tcyt';
         112,111, 1111,   1112,  1113};
    LymphoNets.thresh.FigFlag = 1; %output scatterplot of the cells and nets


    LymphoNets.Size      = uint16(zeros(size(MorpResults.Indexes)));
    LymphoNets.DegreCent = uint16(zeros(size(MorpResults.Indexes)));
    LymphoNets.ClnesCent = uint16(zeros(size(MorpResults.Indexes)));
    LymphoNets.NetworkID = uint16(zeros(size(MorpResults.Indexes)));
    max_inter = 15;
    Interactors = uint32(zeros(size(MorpResults.Indexes,1),max_inter))+NaN;

    summaryNames = {'Tissue ID','Cell #',LymphoNets.thresh.CellTypes{1,:},'Median X','Median Y'};
    LymphoNets.Summary = array2table(zeros(size(summaryNames)),'VariableNames',summaryNames); 

    % make the cell type filter 
    CTfilter = zeros(size(MorpResults.Indexes));
    for f = 1:size(LymphoNets.thresh.CellTypes,2)
        CT = LymphoNets.thresh.CellTypes{2,f};
        layer = ceil(log10(CT+1));
        CTfilter = CTfilter | CellType.Matrix(:,layer) == CT;
    end

    % loop around each tissue
    for j = 1:length(options.MouseNum)    
        disp(['Processing networks in tissue ' filename.tissues{j}])

        % find cells from single tissue and their x-y coordinates
        index = MorpResults.Indexes == j & ROIResults.SpleenIndex == 0;
        coords = double([MorpResults.X(index) MorpResults.Y(index)]);

        % perform delaunay graph for all the cells
        S = delaunaygraph(coords,LymphoNets.thresh.dist,0);

        S_g = graph(S);
    %     figure
    %     plot(S_g,'XData',MorpResults.X(index),'YData',MorpResults.Y(index))

        conn_mat = zeros(sum(index),max_inter);
        for node = 1:size(S_g.Nodes,1)
            temp_neigh = neighbors(S_g,node);
            conn_mat(node,1:length(temp_neigh)) = temp_neigh';
        end
        Interactors(index,:) = conn_mat;
        clear S_g

        % filter the nodes/edges on desired celltypes from CTfilter
        celltypefilt = CellType.Matrix(index,:);
        index_t = CTfilter(index);

        S(~index_t,:) = 0;
        S(:,~index_t) = 0;
        T = graph(S);

        [bins,binsizes]=conncomp(T);
        nets = find(binsizes>LymphoNets.thresh.netsize);

        % plot Lymphonets in space if desired
        if LymphoNets.thresh.FigFlag
            figure
            scatter(double(MorpResults.X(index)),-double(MorpResults.Y(index)),5,'filled','MarkerFaceColor',[.9,.9,.9])
            hold on
            scatter(double(MorpResults.X(index&CTfilter)),-double(MorpResults.Y(index&CTfilter)),5,'filled','MarkerFaceColor',[0.7,0.6,0.85])

            binFilt = ismember(bins,nets);
            plot(subgraph(T,binFilt),'XData',coords(binFilt,1),'YData',-coords(binFilt,2)) % can change the color if desired, but i like that orange ¯\_(?)_/¯
            title(filename.tissues{j})
            legend('all cells','cells filtered on','networks')
        end

        % initialize variables
        degree_cent = zeros(size(bins));
        closen_cent = zeros(size(bins));
        netsize =  zeros(size(bins));
        netId = zeros(size(bins));
        TempSummary = zeros(length(nets),length(summaryNames));

        for n = 1:length(nets)
            % find nodes of subgraph
            S_sub = subgraph(T,bins==nets(n));
            degree_cent(bins==nets(n)) = centrality(S_sub,'degree');
            closen_cent(bins==nets(n)) = centrality(S_sub,'closeness');
            netsize(bins==nets(n)) = binsizes(nets(n));
            netId(bins==nets(n)) = n;
            x_temp = coords(bins==nets(n),1);
            y_temp = coords(bins==nets(n),2);

            TempSummary(n,1) = j;
            TempSummary(n,2) = binsizes(nets(n));
            for CT = 1:size(LymphoNets.thresh.CellTypes,2)
                CTcode = LymphoNets.thresh.CellTypes{2,CT};
                CTlayer = ceil(log10(CTcode+1));
                TempSummary(n,CT+2) = sum(celltypefilt(bins==nets(n),CTlayer) == CTcode);
            end
            TempSummary(n,size(LymphoNets.thresh.CellTypes,2)+3) = median(x_temp);
            TempSummary(n,size(LymphoNets.thresh.CellTypes,2)+4) = median(y_temp);
        end

        LymphoNets.Summary = [LymphoNets.Summary; 
            array2table(TempSummary,'VariableNames',summaryNames)];

        LymphoNets.Size(index) = padarray(netsize,[0 sum(index)-length(netsize)],'post');
        LymphoNets.DegreCent(index) = padarray(degree_cent,[0 sum(index)-length(degree_cent)],'post');
        LymphoNets.ClnesCent(index) = padarray(closen_cent,[0 sum(index)-length(closen_cent)],'post');
        LymphoNets.NetworkID(index) = padarray(netId,[0 sum(index)-length(netId)],'post');

    end
    LymphoNets.Summary = LymphoNets.Summary(2:end,:);

    disp('LymphoNet exploration done')
    %%% SAVE
    save([filename.analfolder filename.resufolder 'Results_Nets_' options.date '_dist' num2str(LymphoNets.thresh.dist) '.mat'],'LymphoNets','Interactors')
    disp(['Done d ' num2str(d) ' of ' num2str(length(distance_cutoff))])
end

%%
thresh = [];
thresh.dist = distance_cutoff(2) ;
thresh.netsize = 5;
ComboCT = { [112, 111, 1111,   1112,  1113 ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [112,      1111,   1112,  1113 ] ; 
            [          1111,   1112,  1113 ]; 
            [112,              1112,  1113 ]
            [112,      1111,          1113 ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [112,      1111,   1112,       ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [112,      1111,   1112,  1113, 121 ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [112,      1111,   1112,  1113, 122 ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [112,      1111,   1112,  1113, 123 ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [112,      1111,   1112,  1113, 124 ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [112,      1111,   1112,  1113, 125 ] ; % 'B','T','Treg','Thelp','Tcyt'; 
            [121, 122, 123, 124, 125 ]
             };
         


for j = 1:length(options.MouseNum)    
    disp(['Processing networks in tissue ' filename.tissues{j}])

    
    
    % find cells from single tissue and their x-y coordinates
    index = MorpResults.Indexes == j & ROIResults.SpleenIndex == 0;
    coords = double([MorpResults.X(index) MorpResults.Y(index)]);

    % perform delaunay graph for all the cells
    S = delaunaygraph(coords,thresh.dist,0);

    
    
    
    combo_names = {};
    % filter the nodes/edges on desired celltypes from CTfilter
    celltypefilt = CellType.Matrix(index,:);
    for combo = 1:length(ComboCT)
        S_type = S;
        
        
        CTfilter = zeros(length(index),1);
        for f = 1:size(ComboCT{combo},2)
            CT = ComboCT{combo}(f);
            layer = ceil(log10(CT+1));
            CTfilter = CTfilter | CellType.Matrix(:,layer) == CT;
            combo_names{combo}{f} = CellType.names{find(CellType.codes == ComboCT{combo}(f))};
        end
        
        if j == 1
            summaryNames{combo} = [{'Tissue ID'},{'Cell #'}, combo_names{combo},{'Median X'},{'Median Y'}];
            Summary{combo}     = array2table(zeros(size(summaryNames{combo})),'VariableNames',summaryNames{combo}); 
            Size{combo} = uint16(zeros(size(MorpResults.Indexes)));
        end
        
        index_t = CTfilter(index);

        S_type(~index_t,:) = 0;
        S_type(:,~index_t) = 0;
        T = graph(S_type);

        [bins,binsizes]=conncomp(T);
        nets = find(binsizes>thresh.netsize);
        
        netsize =  zeros(size(bins));
        TempSummary = zeros(length(nets),length(summaryNames{combo}));
        
        for n = 1:length(nets)
            % find nodes of subgraph
            S_sub = subgraph(T,bins==nets(n));
            netsize(bins==nets(n)) = binsizes(nets(n));
            netId(bins==nets(n)) = n;
            x_temp = coords(bins==nets(n),1);
            y_temp = coords(bins==nets(n),2);

            TempSummary(n,1) = j;
            TempSummary(n,2) = binsizes(nets(n));
            for f = 1:length(ComboCT{combo})
                CTcode = ComboCT{combo}(f);
                CTlayer = ceil(log10(CTcode+1));
                TempSummary(n,f+2) = sum(celltypefilt(bins==nets(n),CTlayer) == CTcode);
            end
            TempSummary(n,length(ComboCT{combo})+3) = median(x_temp);
            TempSummary(n,length(ComboCT{combo})+4) = median(y_temp);
        end

        Summary{combo} = [Summary{combo}; array2table(TempSummary,'VariableNames',summaryNames{combo})];

        Size{combo}(index) = padarray(netsize,[0 sum(index)-length(netsize)],'post');
    end
end

tissues = {1:6,14:21, [1:6 14:21] };

for s = 1:6%length(Summary)
    for t = 1:3%length(tissues{s})
        ind_net = ismember(Summary{1,s}.("Tissue ID"),tissues{t});
        cumnets(s,t) = sum(ind_net);
        meansize(s,t) = mean(Summary{1,s}.("Cell #")(ind_net));
        prc25(s,t) = prctile(Summary{1,s}.("Cell #")(ind_net),25);
        prc75(s,t) = prctile(Summary{1,s}.("Cell #")(ind_net),75);
        
        figure(t)
        subplot(2,3,s)
        [n,h]=hist(Summary{1,s}.("Cell #")(ind_net),linspace(5,40,16)+1.25);
        bar(h,n)
        ylim([0 600])
        sum(Summary{1,s}.("Cell #")(ind_net) <7.5)
    end
end
