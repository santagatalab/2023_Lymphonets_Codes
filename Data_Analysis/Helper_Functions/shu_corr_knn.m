function [corrfunc,rad_approx] = shu_corr_knn(xydata,magvalA,magvalB,k_val)    
    
    magvalA = double(magvalA);
    magvalB = double(magvalB);
    
    % mean center and stdev norm the variables
    magvalA = (magvalA-mean(magvalA))./std(magvalA);
    magvalB = (magvalB-mean(magvalB))./std(magvalB);
    
  
    [idx,dist] = knnsearch(xydata,xydata,'k',k_val,'NSMethod','kdtree');
    rad_approx = mean(dist); %i.e. the k'th neighbor is on average this distance
%     size(magvalA(idx))
%     size(repmat(magvalB,1,size(idx,2)))
    magdots = magvalA(idx).*repmat(magvalB,1,size(idx,2));
    corrfunc = mean(magdots);
end