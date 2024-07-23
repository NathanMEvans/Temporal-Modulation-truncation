%% file locations
data_location = 'Anonymised_Data/imprints';

%% add to the path
addpath("truncations")
addpath("BrewerMap")
%% get patient lists
patients_dir = dir(data_location);
patients = string({patients_dir(3:end).name});


%% create output table
truncation_pairs = table('Size',[0,2],'VariableNames',{'patient_ID','pair_indexes'},'VariableTypes',{'string','cell'});
metadata_table = table();

%% loop over patients
for pat = 1:length(patients)
    %% load data
    patient = patients{pat};
    patient = patient(1:end-4);
 
    load(sprintf("Anonymised_Data/imprints/%s",patient),'imprints','data_tbl')
    if length(imprints) < 2
        continue;
    end

    pairs = find_truncated_pairs(imprints);
    if isempty(pairs)
        continue
    end
    
    % store pairs
    truncation_pairs.patient_ID(end+1) = patient;
    truncation_pairs.pair_indexes(end) = {pairs};

    if pat == 1
        metadata_table = data_tbl;
    else
        metadata_table = cat(1,metadata_table,data_tbl);
    end

end
%% add medication load to seizure metadata
metadata_table.rounded_start = dateshift(metadata_table.start,'start','minute');
metadata_table.combined_drug_load = NaN(height(metadata_table),1);

for pat = 1:height(truncation_pairs)
    
    patient = truncation_pairs.patient_ID{pat};
    figure(pat);
    sgtitle(patient);
    % load BPCs
    csv = sprintf("Anonymised_Data/drug_plasma_concentrations/%s.mat",patient);
    if ~exist(csv,"file")
        continue;
    end
    load(csv,'drug_plasma_concentrations');
    
    drug_names = drug_plasma_concentrations.Properties.VariableNames;
    drug_names = drug_names(2:end); % exclude datetime
    relative_DPC = drug_plasma_concentrations;
    for drug = 1:length(drug_names)
        
        relative_DPC.(drug_names{drug}) = relative_DPC.(drug_names{drug})./(max(relative_DPC.(drug_names{drug})));
        relative_DPC.(drug_names{drug}) = movmean(relative_DPC.(drug_names{drug}),1440);% 12 hour smoothing between max and zero
    end
    relative_DPC.Combined = mean(relative_DPC(:,2:end),2).mean;
    metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)) = relative_DPC.Combined(find(ismember(relative_DPC.DateTime, metadata_table.rounded_start(strcmp(metadata_table.patient_id,patient)))));


    subplot(3,1,1)
    hold on
    plot(relative_DPC,'DateTime','Combined')
    scatter(metadata_table.rounded_start(strcmp(metadata_table.patient_id,patient)),metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)),'filled','red')
    ylim([0,1])
    title('Drug load at seizure events')
end

save("truncation_pairs_metadata","metadata_table","truncation_pairs")
%% Figure showing difference between each pair
load("truncation_pairs_metadata","metadata_table","truncation_pairs")

load_diffs = [];
patients = unique(metadata_table.patient_id);
for pat = 1:length(patients)
    mask = strcmp(metadata_table.patient_id,patients(pat));
    
    pat_drug_load = metadata_table.combined_drug_load(mask);
    figure(pat);
    subplot(3,2,3)
    drug_load_differences = pat_drug_load - pat_drug_load';
    imagesc(drug_load_differences);
    set(gca, 'YDir','normal')
    title("Difference in ASM loads")
    xlabel("Seizures")
    ylabel("Seizures")
    colorbar
    cmap = brewermap(256,'RdBu');
    colormap(cmap);

    pairs =  truncation_pairs.pair_indexes{strcmp(string(patients(pat)),truncation_pairs.patient_ID)};
    pat_load_diffs = [];
    pat_trunc = [];
    pat_cont = [];
    for pair = 1:size(pairs,1)
        continuing = pairs(pair,1);
        truncated = pairs(pair,2);
        rectangle('Position',[continuing-0.5 truncated-0.5 1 1])
        pat_load_diffs(end+1) = pat_drug_load(truncated) -  pat_drug_load(continuing);
    end
    load_diffs(end+1) = median(pat_load_diffs);

    subplot(3,2,4)
    boxplot(pat_load_diffs);
    hold on
    colour_indexes = ceil((pat_load_diffs - min(drug_load_differences,[],"all"))/range(drug_load_differences,"all")*height(cmap));
    scatter(ones(size(pat_load_diffs)),pat_load_diffs,80,cmap(colour_indexes,:),'filled','square','MarkerEdgeColor','b')
    title("truncated pairs")
    ylabel("Difference in ASM loads")
end
%%
for pat = 1:length(patients)
    figure(pat)
    subplot(3,1,3)
    hist = histogram(load_diffs,'BinEdges',-1:0.1:1);
    title("Histogram of difference in medication load between trucnated and continuing seizures for each patient")
    xlabel("medication load difference")
    ylabel("Patient count")
    xline(mean(load_diffs),'red')
    [H,P,CI,STATS] = ttest(load_diffs);
    text(median(load_diffs),max(hist.Values) + 0.5,sprintf("P=%0.2f T=%0.2f",P,STATS.tstat))
    ylim([0,max(hist.Values) + 1])
end

