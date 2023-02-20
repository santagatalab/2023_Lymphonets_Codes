function [all_bins] = FindBins(reference, interest, centroids, max_r, dr)
    all_bins = zeros(length(reference), ceil(max_r/dr));
    for j=1:length(reference)
        curr_centroid = reference(j,:);
        crop_idcs = find(abs(centroids(:,1)-curr_centroid(1)) < max_r & abs(centroids(:,2)-curr_centroid(2)) < max_r);
        cropped_centroids = centroids(crop_idcs,:);

        idcs = find(abs(interest(:,1)-curr_centroid(1)) < max_r & abs(interest(:,2)-curr_centroid(2)) < max_r);
        filtered = interest(idcs,:);
        all_bins(j,:) = RadialDistribution(curr_centroid, filtered, cropped_centroids, max_r, dr);

%         figure
%         scatter(cropped_centroids(:,1), cropped_centroids(:,2), 3, 'filled')
%         hold on
%         scatter(filtered(:,1), filtered(:,2), 8, 'm', 'filled')
%         hold on
%         scatter(curr_centroid(:,1), curr_centroid(:,2), 50, 'g', 'x')
%         hold on
%         for k=1:ceil(max_r/dr)
%             viscircles(curr_centroid, dr*k, 'LineWidth', 1, 'LineStyle', '--');
%             hold on
%         end
%         axis equal
%         keyboard
    end
end