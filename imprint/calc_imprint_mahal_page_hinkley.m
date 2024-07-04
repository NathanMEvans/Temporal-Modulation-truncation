function [data_tbl, cell_imprint] = calc_imprint_mahal_page_hinkley(data_tbl, opts)%
% Compute seizure imprint based on EEG recordings

% input:
%   - data_tbl: full data table
%   - optional inputs
%       - window_size: window size for which imprint will be computed
%       - min_sz_count: minimum number of seizures to have been recorded per patient
%       - folder: folder to store markers in

% output
%   - data_tbl: full data table
%   - cell_imprint

    arguments
        data_tbl
        opts.window_size (1,1) double {mustBeNumeric} = 1; % Decide window size (in seconds) for which markers and imprint are computed
        opts.folder = 'onset_calcs'; % folder to store markers in
        opts.window_overlap (1,1) double {mustBeNumeric} = 7/8; % amount that windows overlap (in seconds) - must be less than window size
        opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5
        opts.pre_buffer (1,1) double {mustBeInteger, mustBePositive} = 10 %in seconds
        opts.ict_buffer (1,1) double {mustBeNumeric} = 10 % Number of seconds to shift the clinicalaly marked onset 
        opts.lambda (1,1) double {mustBeNumeric, mustBePositive} = 200 %
    end
    
    %fill in optional arguments
    window_size = opts.window_size;
    folder = opts.folder;
    window_overlap = opts.window_overlap;
    mad_thresh = opts.mad_thresh;
    pre_buffer = opts.pre_buffer;
    ict_buffer = opts.ict_buffer;
    lambda = opts.lambda;
    
    % Set basefolder to store markers
    patient = data_tbl.patient_id{1};
    basefolder = sprintf('%s/%s', folder, string(patient));
    sampling_rate=data_tbl.segment_fs(1) ;%just using the first, as all sampling was the same in our data after preproc.
    
    linelength_db = LL_db([basefolder '/LL_db/']); %setup folder for all Line Length measures
    linelength_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    linelength_db.paramset_tbl                       % display all currently tracked paramsets
    linelength_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    energy_db = Energy_db([basefolder '/Energy_db']);
    energy_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    energy_db.paramset_tbl                       % display all currently tracked paramsets
    energy_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    %bands: [1 4; 4 8; 8 13; 13 30; 30 60, 60 100]; 
    bandpower_db = BP_db([basefolder '/BP_db']);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[1 4]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[4 8]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[8 13]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[13 30]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[30 60]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[60 100]);
    bandpower_db.paramset_tbl                       % display all currently tracked paramsets
    bandpower_db.calc(data_tbl,[]);   
    
    
    %% Import markers (line length, energy, bandpowers)
    calcs_ll = linelength_db.get(1);    % get all calculation outputs in variable. 
    calcs_energy = energy_db.get(1);
    calcs_bp_delta = bandpower_db.get(1);
    calcs_bp_theta = bandpower_db.get(2);
    calcs_bp_alpha = bandpower_db.get(3);
    calcs_bp_beta = bandpower_db.get(4);
    calcs_bp_gamma = bandpower_db.get(5);
    calcs_bp_hgamma = bandpower_db.get(6);
    
    %% calculate imprint and recruitment markers
    val_tbl=[calcs_ll.LL_ms calcs_energy.energy_ms calcs_bp_delta.bp ...
        calcs_bp_theta.bp calcs_bp_alpha.bp calcs_bp_beta.bp calcs_bp_gamma.bp calcs_bp_hgamma.bp];
    %% Create a table of seizure data
    sz_mat_tab = table(data_tbl.segment_id, cell(size(val_tbl,1),1), calcs_ll.t_wndw,...
        'VariableNames',["segment_id", "feat_mat", "tw"]);
    for sz = 1:size(val_tbl,1)
        sz_mat=log(cat(3,val_tbl{sz,:})); % Here we log-transform markers 
        sz_mat_tab.feat_mat{sz} = sz_mat;
    end
    [imprint_out,cell_imprint,~,cell_madscores, cell_pre_features_mad,...
        cell_pre_mahal_mat] = mahal_imprint_page_hinkley(...
            data_tbl, ...
            sz_mat_tab,...
            calcs_ll.t_wndw, ...
            'mad_thresh',mad_thresh,...
            "lambda",lambda, ...
            "pre_buffer", pre_buffer,...
            "ict_buffer", ict_buffer...
         );  % using the same t_wndw for all features as using same window length or overlap
    
    cell_imprint.cell_madscores = cell_madscores;  
    cell_imprint.cell_pre_features_mad = cell_pre_features_mad;
    cell_imprint.cell_pre_mahal_mat = cell_pre_mahal_mat;
    cell_imprint.onset_time = imprint_out.("Onset time") - ict_buffer;
end