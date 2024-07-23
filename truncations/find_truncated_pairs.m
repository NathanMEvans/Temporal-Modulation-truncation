function [pairs, ratio_change] = find_truncated_pairs(imprints)
    %% cut at intial zero and maximal spread
    onset_t = zeros(size(imprints));
    max_spread_t = zeros(size(imprints));
    max_spread_n = zeros(size(imprints));
    for i = 1:length(imprints)
        imprint = imprints{i};
        onset_t(i) = find(any(imprint,1),1);
        [max_spread_n(i), max_spread_t(i)] = max(sum(imprint,1,'omitmissing'));
        imprints{i} = imprint(:,onset_t(i):max_spread_t(i));
    end
    %% then find difference with DTW
    [dtw_oe_dist_norm,dtw_dist_norm, fraction_oe, ~] = all_dtw_oe_distances(imprints);
    ratio_change = ((dtw_dist_norm - dtw_oe_dist_norm) ./ (dtw_dist_norm));

    %%
    ratio_change(fraction_oe < 10) = 0; % ignore if less than 10% of seizure is used
    ratio_change(ratio_change < 0.2) = 0; % ignore if less than 20% impoved on non truncated comparison
    ratio_change(max_spread_n*(max_spread_n.^-1)' < 1.5) = 0; % ignore if seizure doesn't spread to 50% more channels than truncated

    if ~any(ratio_change,"all")
        pairs = [];
        return;
    end
    [continuing_seizures,truncated_seizures] = find(ratio_change);
    pairs = [continuing_seizures,truncated_seizures];
end