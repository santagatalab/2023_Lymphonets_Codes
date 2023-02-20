function [bins] = RadialDistribution(curr_centroid, filtered_centroids, centroids, max_r, dr)
    D_filtered = pdist2(curr_centroid, filtered_centroids, 'euclidean');
    D_filtered = D_filtered(D_filtered>0 & D_filtered<max_r);

    D_all = pdist2(curr_centroid, centroids, 'euclidean');
    D_all = D_all(D_all>0 & D_all<max_r);

    num_bins = ceil(max_r/dr);
    bins = zeros(num_bins,1);
    for i=1:num_bins
        cutoff_lo = (i-1)*dr;
        cutoff_hi = i*dr;
        filtered_count = sum(D_filtered>=cutoff_lo & D_filtered<cutoff_hi);
        tot_count = sum(D_all>=cutoff_lo & D_all<cutoff_hi);
        if tot_count == 0
            bins(i) = 0;
        else
            bins(i) = filtered_count/tot_count;
        end
    end
end