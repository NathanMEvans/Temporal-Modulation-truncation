function [dtw_oe_dist_norm,dtw_dist_norm, fraction_oe, truncation_points] = all_dtw_oe_distances(data)
    dtw_oe_dist = zeros(length(data));
    fraction_oe = zeros(length(data));
    truncation_points = zeros(length(data));
    dtw_dist = zeros(length(data));
    dtw_oe_dist_norm = zeros(length(data));
    dtw_dist_norm = zeros(length(data));
    for i=1:length(data)
        data_i = data{i}; 
        parfor j = 1:length(data) % could be parfor
            [dtw_oe_dist(i,j), fraction_oe(i,j), truncation_points(i,j), dtw_dist(i,j)]=dtw_openEnded(data_i',data{j}');
            dtw_oe_dist_norm(i,j)=dtw_oe_dist(i,j)./max(size(data_i,2),size(data{j},2));
            dtw_dist_norm(i,j)=dtw_dist(i,j)./max(size(data_i,2),size(data{j},2));
        end
    end
end

