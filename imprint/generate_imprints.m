%% script to generate imprints from EEG data
% EEG data should be saved in the Data/EEG_dataexport directory, or the
% value for the data_location variable should be changed
% EEG data should be saved in files names [patientID].mat
% the files should have a single data_export table variable containing one
% row for each recorded seizure and requieres the following columns
%   segment_id      a unique identifier for the recording
%   patient_id      a unique identifier for the patient
%   duration        the duration of the recorded seizure in seconds
%   segment_pre     the length of the recording included before seizure onset in seconds - recomended to be at 120 seconds
%   segment_post    the length of the recording included after seizure conclusion in seconds
%   segment_fs      the frequency of the recording in hz
%   segment_data    the recording data in a channels by time points matrix

date_format = 'yyyy-MM-dd''T''HH:mm:ss''Z''';

%% file locations
data_location = '../Data/EEG_dataexport/';


%% add to the path
addpath(genpath('../ictal_onset/')); % calc_imprint_mahal_page_hinkley

%% get patient lists
patients_dir = dir(data_location);
patients = string({patients_dir(3:end).name});

%% create output table
imprints = table('Size',[0,2],'VariableNames',{'patient_ID','imprints'},'VariableTypes',{'string','cell'});
metadata_table = table();

%% loop over patients
for pat = 1:length(patients)
    %% load data
    patient = patients{pat};
    patient = patient(1:end-4);
    data_export = {};
    load(sprintf('%s/%s.mat', data_location, patient), 'data_export');
    if isempty(data_export)
        continue
    end
    %% calculate imprint
    mad_thresh = 5;
    rec_thresh = 32;
    lambda = 200;
    prop_rec = 0.7;
    [data_tbl, cell_imprint_mahal] = calc_imprint_mahal_page_hinkley(data_export,...
        'window_size', 1, ...
        'window_overlap', 7/8, ...
        'folder','tmp', ...%'folder','../data/onset_calcs', ...
        "mad_thresh",mad_thresh, ...
        "lambda",lambda...
        );
    imprints = cell_imprint_mahal.cell_imprint;
    
     %% remove empty imprints
    remove = false(size(imprints));
    for i = 1:length(imprints)
        imprint = imprints{i};
        if ~any(imprint,"all")
            remove(i) = true;
        end
    end
    imprints(remove) = [];
    data_tbl(remove,:) = [];
    data_tbl.segment_data = [];
    save("../Data/imprints/"+patient,"data_tbl","imprints");
end