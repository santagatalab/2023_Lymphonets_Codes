function [kvals, lambda] = ripleys(xydata, rvals, subsample)    

    if subsample && (length(xydata) > 50000) % take bounding box that is 1/4 the size of the old one
        xmin = (min(xydata(:,1)) + max(xydata(:,1)))*(1/4);
        xmax = (min(xydata(:,1)) + max(xydata(:,1)))*(3/4);
        ymin = (min(xydata(:,2)) + max(xydata(:,2)))*(1/4);
        ymax = (min(xydata(:,2)) + max(xydata(:,2)))*(3/4);
        xfilt = xydata(:,1) > xmin & xydata(:,1) < xmax;
        yfilt = xydata(:,2) > ymin & xydata(:,2) < ymax;
        xydata = xydata(xfilt & yfilt, :);      
        disp(['new cell count: ' num2str(length(xydata))])
        
        if length(xydata) > 50000
            xmin = (min(xydata(:,1)) + max(xydata(:,1)))*(1/4);
            xmax = (min(xydata(:,1)) + max(xydata(:,1)))*(3/4);
            ymin = (min(xydata(:,2)) + max(xydata(:,2)))*(1/4);
            ymax = (min(xydata(:,2)) + max(xydata(:,2)))*(3/4);
            xfilt = xydata(:,1) > xmin & xydata(:,1) < xmax;
            yfilt = xydata(:,2) > ymin & xydata(:,2) < ymax;
            xydata = xydata(xfilt & yfilt, :);      
            disp(['new cell count 2: ' num2str(length(xydata))])
        end
    end
    kvals = [];
    [~,D] = rangesearch(xydata,xydata,max(rvals),'NSMethod','kdtree','SortIndices',false);
    D = [D{:}];
    n = double(size(xydata,1));
    lambda = n/(double(range(xydata(:,1)))*double(range(xydata(:,2))));
    
    for r=rvals        
        K = (sum(D(D~=0)<r)/n)*(1/lambda);
        kvals = [kvals K];
    end
end






