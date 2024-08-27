%% file locations
data_location = 'Anonymised_Data/imprints';

%% add to the path
addpath("truncations")
addpath("BrewerMap")
%% get patient lists
patients_dir = dir(data_location);
patients = string({patients_dir(3:end).name});


%% create output table
truncation_pairs = table('Size',[0,6],'VariableNames',{'patient_ID','pair_indexes','n_seizures','n_truncated_seizures','n_continuing_seizures','n_involved_seizures'},'VariableTypes',{'string','cell','double','double','double','double'});
metadata_table = table();

%% loop over patients
for pat = 1:length(patients)
    fprintf("Subject %d \n", pat)
    %% load data
    patient = patients{pat};
    patient = patient(1:end-4);
 
    load(sprintf("Anonymised_Data/imprints/%s",patient),'imprints','data_tbl')
    
    truncation_pairs.patient_ID(end+1) = patient;
    if length(imprints) > 1
        % store pairs
        pairs = find_truncated_pairs(imprints);
        truncation_pairs.pair_indexes(end) = {pairs};

        % count truncations
        if ~isempty(pairs)
            truncation_pairs.n_seizures(end) = length(imprints);
            truncation_pairs.n_truncated_seizures(end) = length(unique(pairs(:,2)));
            truncation_pairs.n_continuing_seizures(end) = length(unique(pairs(:,1)));
            truncation_pairs.n_involved_seizures(end) = length(unique([pairs(:,1),pairs(:,2)]));
        end
    end
   
    if size(metadata_table,1) == 0
        metadata_table = data_tbl;
    else
        metadata_table = cat(1,metadata_table,data_tbl);
    end

end
%% add medication load to seizure metadata
metadata_table.rounded_start = dateshift(metadata_table.start,'start','minute');
metadata_table.combined_drug_load = NaN(height(metadata_table),1);
has_tapering = [];
for pat = [1:height(truncation_pairs)]
    
    patient = truncation_pairs.patient_ID{pat};
    figure(pat);
    sgtitle(patient);
    % load BPCs
    csv = sprintf("Anonymised_Data/drug_plasma_concentrations/%s.mat",patient);
    if ~exist(csv,"file")
        continue;
    end
    load(csv,'drug_plasma_concentrations');
    has_tapering(end+1) = pat;
    drug_names = drug_plasma_concentrations.Properties.VariableNames;
    drug_names = drug_names(2:end); % exclude datetime
    relative_DPC = drug_plasma_concentrations;
    for drug = 1:length(drug_names)
        
        mov_mean = movmean(relative_DPC.(drug_names{drug}),1440*2); % 12*2 hour smoothing (30 seconds per row in CSV)
        mov_mean([1:720*2, end-720*2:end]) = NaN;
        relative_DPC.(drug_names{drug}) = mov_mean;
        relative_DPC.(drug_names{drug}) = relative_DPC.(drug_names{drug})./max(relative_DPC.(drug_names{drug}),[],'omitnan'); % normalise max to 1
    end
    relative_DPC.Combined = mean(table2array(relative_DPC(:,2:end)),2);
    steady_state = mean(relative_DPC.Combined(721*2:5760));
    relative_DPC.Combined = relative_DPC.Combined./steady_state;
    metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)) = relative_DPC.Combined(find(ismember(relative_DPC.DateTime, metadata_table.rounded_start(strcmp(metadata_table.patient_id,patient)))));

    subplot(10,10,[1:8,11:18])
    hold on
    plot(relative_DPC,'DateTime','Combined')
    scatter(metadata_table.rounded_start(strcmp(metadata_table.patient_id,patient)),metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)),'filled','black')
    ylim([0,max(relative_DPC.Combined)])
    %title('Drug load at seizure events')
    xlabel("");
    ax = gca;
    ax.XTickLabel = ax.XTickLabel;

    n_seizures = height(metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)));
    subplot(10,10,[22:25])
    plot(metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)))
    hold on
    scatter(1:n_seizures,metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)),'filled','black')
    xlim([0.5,n_seizures+0.5])
    xticks([])

    subplot(10,10,[31,41,51,61])
    plot(metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)))
    hold on
    scatter(1:n_seizures,metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)),'filled','black')
    xlim([0.5,n_seizures+0.5])
    view([90 -90])
    set(gca, 'xDir','reverse')
    set(gca, 'YDir','reverse')
end

save("truncation_pairs_metadata","metadata_table","truncation_pairs")

%% Figure showing difference between each pair
load("truncation_pairs_metadata","metadata_table","truncation_pairs")

patients = unique(metadata_table.patient_id);
load_diffs = nan(1,length(patients));

for pat = [1:height(truncation_pairs)]
    mask = strcmp(metadata_table.patient_id,patients(pat));
    
    pat_drug_load = metadata_table.combined_drug_load(mask);
    if all(isnan(pat_drug_load))
        continue
    end
    figure(pat);
    subplot(10,10,[32:36,42:46,52:56,62:66])
    drug_load_differences = pat_drug_load - pat_drug_load';
    imagesc(drug_load_differences);
    set(gca, 'YDir','reverse')
    yticks([])
    %title("Difference in ASM loads")
    xlabel("Seizures")
    %ylabel("Seizures")
    colorbar
    cmap = brewermap(256,'RdBu');
    colormap(cmap);

    pairs =  truncation_pairs.pair_indexes{strcmp(string(patients(pat)),truncation_pairs.patient_ID)};
    if isempty(pairs)
        continue;
    end
    pat_load_diffs = [];
    pat_trunc = [];
    pat_cont = [];
    for pair = 1:size(pairs,1)
        continuing = pairs(pair,1);
        truncated = pairs(pair,2);
        rectangle('Position',[continuing-0.5 truncated-0.5 1 1],'LineWidth',2)
        pat_load_diffs(end+1) = pat_drug_load(truncated) -  pat_drug_load(continuing);
    end
    load_diffs(pat) = median(pat_load_diffs);


    subplot(10,10,[28,38,48,58,68])
    
    boxchart(ones(size(pat_load_diffs)),pat_load_diffs);
    hold on
    colour_indexes = ceil((pat_load_diffs - min(drug_load_differences,[],"all"))/range(drug_load_differences,"all")*height(cmap));
    scatter(ones(size(pat_load_diffs)),pat_load_diffs,80,cmap(colour_indexes,:),'filled','square','MarkerEdgeColor','b')
    
    %title("truncated pairs")
    ylabel("Difference in ASM loads")
    ylim([-1,1])
end
%%
for pat = [1:height(truncation_pairs)]
    figure(pat)
    subplot(10,10,[71:78,81:88,91:98])


    hist = histogram(load_diffs,'BinEdges',-1:0.1:1);
    hist_values = hist.Values;
    xs = -0.95:0.1:0.95;
    for n_x = 1:20
        for y = 1:hist_values(n_x)
            scatter(xs(n_x),y,'Black','filled')
            hold on
        end
    end
    hold off
    xlim([-1.01,1.05])

    %title("Histogram of difference in medication load between truncated and continuing seizures for each patient")
    xlabel("medication load difference")
    ylabel("Patient count")
    xline(mean(load_diffs, "omitnan"),'red')
    [H,P,CI,STATS] = ttest(load_diffs);
    text(median(load_diffs, 'omitnan'),max(hist_values) + 0.5,sprintf("P=%0.2f T=%0.2f",P,STATS.tstat))
    ylim([0.5,max(hist_values) + 1])
end

fprintf("mean = %f\n",mean(load_diffs,'omitmissing'))
fprintf("std = %f\n",std(load_diffs,'omitmissing'))
fprintf("p = %f\n",P)
STATS
% 
% 
% figure()
% 
% subplot(5,5,2:5)
% subplot(5,5,[6,11,16,21])
% subplot(5,5,[7:10,12:15,27:20,22:25])
% imagesc(drug_load_differences);
% set(gca, 'YDir','normal')
% title("Difference in ASM loads")
% xlabel("Seizures")
% ylabel("Seizures")
% colorbar
% cmap = brewermap(256,'RdBu');
% colormap(cmap);
% pairs =  truncation_pairs.pair_indexes{strcmp(string(patients(pat)),truncation_pairs.patient_ID)};
% pat_load_diffs = [];
% pat_trunc = [];
% pat_cont = [];
% for pair = 1:size(pairs,1)
% continuing = pairs(pair,1);
% truncated = pairs(pair,2);
% rectangle('Position',[continuing-0.5 truncated-0.5 1 1])
% pat_load_diffs(end+1) = pat_drug_load(truncated) -  pat_drug_load(continuing);
% end
% load_diffs(end+1) = median(pat_load_diffs);
% title()
% title("")
% sgtitle("Difference in ASM loads")
% ylabel("")
% colorbar
% colorbar("off")
% plot(metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)))
% hold on
% scatter(1:height(metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient))),metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)))
% xlim(0.5,23.5)
% xlim([0.5,23.5])
% plot(metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)))
% hold on
% scatter(1:height(metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient))),metadata_table.combined_drug_load(strcmp(metadata_table.patient_id,patient)))
% xlim([0.5,23.5])
% view([90 -90])
% set(gca, 'xDir','reverse')
% set(gca, 'yDir','reverse')
% set(gca, 'xDir','normal')
% set(gca, 'xDir','reverse')
% set(gca, 'YDir','reverse')

%% metadata

grouped_meta = groupsummary(metadata_table,{'patient_id','age','epilepsy_duration','sex','op_type'},'range','combined_drug_load');

grouped_meta_with_truncs = join(grouped_meta,truncation_pairs,'LeftKeys','patient_id','RightKeys','patient_ID');
grouped_meta_with_truncs.has_tapering = ~isnan(grouped_meta_with_truncs.range_combined_drug_load);
grouped_meta_with_truncs.hass_truncations = ~cellfun('isempty', grouped_meta_with_truncs{:,'pair_indexes'} );


age_tapering =  groupsummary(grouped_meta_with_truncs,{'has_tapering',},{'median','std'},{'age','epilepsy_duration'});
age_tapering_truncs =  groupsummary(grouped_meta_with_truncs,{'has_tapering','hass_truncations'},{'median','std'},{'age','epilepsy_duration'});
sex_tapering = groupsummary(grouped_meta_with_truncs,{'has_tapering','sex'});
sex_tapering_truncs = groupsummary(grouped_meta_with_truncs,{'has_tapering','hass_truncations','sex'});
op_tapering = groupsummary(grouped_meta_with_truncs,{'has_tapering','op_type'});
op_tapering_truncs = groupsummary(grouped_meta_with_truncs,{'has_tapering','hass_truncations','op_type'});

age_all = groupsummary(grouped_meta,{},{'median','std'},{'age','epilepsy_duration'});
sex_all = groupsummary(grouped_meta,{'sex'});
op_all = groupsummary(grouped_meta,{'op_type'});

demographics = table('Size',[5,3],'VariableNames',{'All', 'Tapered', 'Tapered with truncated pairs'},'VariableTypes',{'string','string','string'},'RowNames',{'Count', 'Age', 'Sex', 'Disease duration', 'TLE/ETLE'});
demographics{"Count","All"} =               string(sprintf('%i',age_all.GroupCount));
demographics{"Age","All"} =                 string(sprintf('%i(%.1f)',age_all.median_age,age_all.std_age));
demographics{"Sex","All"} =                 string(sprintf('%i/%i',sex_all.GroupCount(sex_all.sex == 'M'),sex_all.GroupCount(sex_all.sex == 'F')));
demographics{"Disease duration","All"} =    string(sprintf('%i(%.1f)',age_all.median_epilepsy_duration,age_all.std_epilepsy_duration));
demographics{"TLE/ETLE","All"} =            string(sprintf('%i/%i',op_all.GroupCount(strcmp(op_all.op_type,'T Lx')),sum(op_all.GroupCount(~strcmp(op_all.op_type,'T Lx')))));

age_tapering = age_tapering(age_tapering.has_tapering,:);
sex_tapering = sex_tapering(sex_tapering.has_tapering,:);
op_tapering = op_tapering(op_tapering.has_tapering,:);

demographics{"Count","Tapered"} =               string(sprintf('%i',age_tapering.GroupCount));
demographics{"Age","Tapered"} =                 string(sprintf('%i(%.1f)',age_tapering.median_age,age_tapering.std_age));
demographics{"Sex","Tapered"} =                 string(sprintf('%i/%i',sex_tapering.GroupCount(sex_tapering.sex == 'M'),sex_tapering.GroupCount(sex_tapering.sex == 'F')));
demographics{"Disease duration","Tapered"} =    string(sprintf('%i(%.1f)',age_tapering.median_epilepsy_duration,age_tapering.std_epilepsy_duration));
demographics{"TLE/ETLE","Tapered"} =            string(sprintf('%i/%i',op_tapering.GroupCount(strcmp(op_tapering.op_type,'T Lx')),sum(op_tapering.GroupCount(~strcmp(op_tapering.op_type,'T Lx')))));

age_tapering_truncs = age_tapering_truncs(age_tapering_truncs.has_tapering & age_tapering_truncs.hass_truncations ,:);
sex_tapering_truncs = sex_tapering_truncs(sex_tapering_truncs.has_tapering & sex_tapering_truncs.hass_truncations ,:);
op_tapering_truncs = op_tapering_truncs(op_tapering_truncs.has_tapering & op_tapering_truncs.hass_truncations ,:);

demographics{"Count","Tapered with truncated pairs"} =               string(sprintf('%i',age_tapering_truncs.GroupCount));
demographics{"Age","Tapered with truncated pairs"} =                 string(sprintf('%i(%.1f)',age_tapering_truncs.median_age,age_tapering_truncs.std_age));
demographics{"Sex","Tapered with truncated pairs"} =                 string(sprintf('%i/%i',sex_tapering_truncs.GroupCount(sex_tapering_truncs.sex == 'M'),sex_tapering_truncs.GroupCount(sex_tapering_truncs.sex == 'F')));
demographics{"Disease duration","Tapered with truncated pairs"} =    string(sprintf('%i(%.1f)',age_tapering_truncs.median_epilepsy_duration,age_tapering_truncs.std_epilepsy_duration));
demographics{"TLE/ETLE","Tapered with truncated pairs"} =            string(sprintf('%i/%i',op_tapering_truncs.GroupCount(strcmp(op_tapering_truncs.op_type,'T Lx')),sum(op_tapering_truncs.GroupCount(~strcmp(op_tapering_truncs.op_type,'T Lx')))));
demographics

