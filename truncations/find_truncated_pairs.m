function pairs = find_truncated_pairs(imprints)
    %% cut at intial zero and maximal spread
    onset_t = zeros(size(imprints));
    max_spread_t = zeros(size(imprints));
    max_spread_n = zeros(size(imprints));
    for i = 1:length(imprints)
        imprint = imprints{i};
        onset_t(i) = find(any(imprint,1),1);
        [max_spread_n(i), max_spread_t(i)] = max(sum(imprint,1,'omitnan'));
        imprints{i} = imprint(:,onset_t(i):max_spread_t(i));
    end
    %% then find difference with DTW
    [dtw_oe_dist_norm,dtw_dist_norm, fraction_oe, ~] = all_dtw_oe_distances(imprints);
    ratioChange = ((dtw_dist_norm - dtw_oe_dist_norm) ./ (dtw_dist_norm));

    %%
    ratioChange(fraction_oe < 10) = 0; % ignore if less than 10% of seizure is used
    ratioChange(ratioChange < 0.2) = 0; % ignore if less than 20% impoved on non truncated comparison
    ratioChange(max_spread_n*(max_spread_n.^-1)' < 1.5) = 0; % ignore if seizure doesn't spread to 50% more channels than truncated

    if ~any(ratioChange,"all")
        pairs = [];
        return;
    end
    [continuing_seizures,truncated_seizures] = find(ratioChange);
    pairs = [continuing_seizures,truncated_seizures];
end