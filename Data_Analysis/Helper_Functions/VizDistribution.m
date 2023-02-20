function VizDistribution(all_bins, smooth_size, plot_size, color_lim, labels, plot_title)
    figure
    for i=1:length(all_bins)
        mouse_bins = all_bins{i};
        pk_locs = zeros(size(mouse_bins,1),1);
        for j=1:size(mouse_bins,1)
            smoothed_bin = movmean(mouse_bins(j,:),smooth_size);
            max_pk_idcs = find(smoothed_bin==max(smoothed_bin));
            pk_locs(j) = max_pk_idcs(1);
        end
        sorted_bins = sortrows([mouse_bins pk_locs], size(mouse_bins,2)+1);
        subplot(plot_size(1),plot_size(2),i)
        imagesc(sorted_bins(:,1:size(mouse_bins,2)));
        caxis(color_lim)
        title(num2str(labels(i)))
    end
    suptitle(plot_title)
%     set(gcf,'color','w');
end