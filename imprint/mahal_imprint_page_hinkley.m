function  [tbl_imprint_out, cell_imprint, cell_t, cell_madscores,cell_pre_features_mad, cell_pre_mahal_mat]  =...
mahal_imprint_page_hinkley(data_tbl,sz_mat_tab,t,opts)

% calculates "imprint" using mahal distance and Page Hinkley style threshold algorithm.
% 
% inputs:
%   - data_tbl: table with columns for ids,pre,duration,fs, and data for better
%       plots
%   - val_cell: cell array with data matrix chns-by-time for each segment,
%       cell columns are different features, rows are different segments
%   - t: cell array with time arrays for each segment where time
%       arrays contain time points in seconds corresponding to the columns of
%       data matrix of that segment
%   - opts: for more info on parameters look below
%outputs are in order of data_tbl, as we assume that input is in order of
%data_tbl

    arguments
        data_tbl
        sz_mat_tab
        t
        opts.mad_thresh (1,1) double {mustBeNumeric, mustBePositive} = 5 %
        % corresponds to the magnitude of changes that should not raise an alarm.
        opts.lambda (1,1) double {mustBeNumeric, mustBePositive} = 200 % 
        % cummilative amount above mad_thresh that signal needs to be to
        % trigger
        opts.pre_buffer (1,1) double {mustBeNumeric, mustBePositive} = 10 %
        opts.ict_buffer% number of seconds before CLO time to go back
    end
    
    pre_buffer=opts.pre_buffer;
    ict_buffer=opts.ict_buffer;
    mad_thresh=opts.mad_thresh;
    lambda = opts.lambda;

    mc=-1/(sqrt(2)*erfcinv(3/2)); % fixed factor for MAD score calculation
    
    if size(data_tbl,1) ~= size(sz_mat_tab,1) || size(data_tbl,1) ~= numel(t) || ...
        size(sz_mat_tab,1) ~= numel(t)
        error('All inputs must have same amount of rows.')
    end
    
    nsegs = size(data_tbl,1);
    cell_madscores = table(cell(nsegs,1), cell(nsegs,1),...
        'VariableNames', ["mahal_MAD", "Un"]);
    cell_pre_mahal_mat = cell(nsegs,1);
    cell_t = cell(nsegs,1);
    cell_pre_features_mad = cell(nsegs,1);

    patient = sprintf("UCLH%s", string(data_tbl(1,:).patient_id));

    %% Create empty output tables
   tbl_imprint_out = table(repelem(patient,nsegs,1), data_tbl.segment_id,...
        nan(nsegs,1), cell(nsegs,1), nan(nsegs,1),...
        'VariableNames',["patient_id", "segment_id", "Onset time", "Onset",...
        "num_chan_imprint"]);

    cell_imprint = table(data_tbl.segment_id, cell(nsegs,1), 'VariableNames',...
       ["segment_id", "cell_imprint"]);

    %% Iterate through seizures 
    for s = 1:nsegs
        fprintf("\n Seizure %d", s)
        %read out seizure segment & compress all features into one 3D array
       
        tw=t{s};%timing for each seizure, in actual seconds
        pre_ids = find(tw<=-(1*pre_buffer+1));
        pre_vals = sz_mat_tab.feat_mat{s}(:,pre_ids,:);
        ictal_ids = find(tw>-1*ict_buffer& tw<=data_tbl.duration(s));
%         ictal_ids = find(tw>-1*ict_buffer& tw<=data_tbl.duration(s));
        ict_vals = sz_mat_tab.feat_mat{s}(:,ictal_ids,:);
        
        %% For each channel create a distribution of Mahalanobis distances (remove
        % time point and compare distance against all other time points)
        % Create a distribution of Mahalanobis distances (one value per
        % channel, marker, and time point)
        pre_mahal_mat = nan(size(pre_vals,1), size(pre_vals(1,:,1),2));
        for chan = 1:size(pre_vals,1)
            pre_mahal_mat(chan, :) = mahal(squeeze(pre_vals(chan,:,:)), squeeze(pre_vals(chan,:,:)));
        end
        
        %% MAD score ictal Mahalanobis distances against preictal Mahalanobis distances
        ict_mahal_mat = nan(size(ict_vals,1), size(ict_vals(1,:,1),2));
        pre_mahal_mad_mat = nan(size(pre_vals,1), size(pre_vals(1,:,1),2));
        mahal_mad_mat = nan(size(ict_vals,1), size(ict_vals(1,:,1),2));
        
        for chan = 1:size(pre_vals,1) 
            pre_dist = pre_mahal_mat(chan,:);
            % Compute preictal median Mpre_mahal_matahalanobis distance
            pre_med = median(pre_dist, 2, 'omitnan');
            pre_smad=mc*mad_rewrite(pre_dist,1,2);
            % Compute preictal MAD values
            pre_mahal_mad_chan = (pre_dist-pre_med)./pre_smad;
            % Remove any preictal windows with MAD >= MAD threshold
            % (essentially removing preictal noise)
            pre_mahal_mad_chan(pre_mahal_mad_chan >= mad_thresh) = NaN;
            
            ref_vals = squeeze(pre_vals(chan,:,:));
            ref_vals(pre_mahal_mad_chan >= mad_thresh,:) = [];

            pre_mahal_mad_mat(chan,:) = pre_mahal_mad_chan;%score preictal to median & scaled mad
            % Remove preictal outliers from pre_dist before computing MAD
            pre_dist = mahal(ref_vals, ref_vals);
            % Recompute med and smad
            pre_med = median(pre_dist, 1, 'omitnan');
            pre_smad=mc*mad_rewrite(pre_dist',1,2);
            % store updated pre_mahal_mat for export 
            pre_mahal_mat(chan,:) = pre_dist;
            % Compute ictal MAD values (ictal values MAD scored against
            % preictal distribution)
            ict_mahal_mat(chan, :) = mahal(squeeze(ict_vals(chan,:,:)), ref_vals);
            mahal_mad_mat(chan,:) = (ict_mahal_mat(chan,:)-pre_med)./pre_smad; %score ictal to median & scaled mad
        end
        
        %% Add inclusion criteria
        v = mad_thresh;
        v_array = ones(size(mahal_mad_mat)) * v;
        Un = cumsum(mahal_mad_mat - v_array,2);

        working_Un = Un; % so we keep a copy of initial 
        t_seg = size(Un,2);
        
        imprint = false(size(mahal_mad_mat));
        % we create a matrix that goes from 1 in column 1 to max_T for
        % column max_T so we can do a comparison to inflection points later
        increasing_t = repmat(1:size(mahal_mad_mat,2),size(mahal_mad_mat,1),1);
        for N = 1:size(mahal_mad_mat,2)
            %% check if above threshold by Lambda
            [un,ind] = min(working_Un(:,1:N),[],2);
            UnDiff = working_Un(:,N) - un; 
            
            pastThreshold_up = UnDiff > lambda; % channels that are pass the threshold
            working_Un(increasing_t < ind & pastThreshold_up) = NaN; % delete points before inflection so downwards threshold can be checked
            
            %% check if bellow threshold by Lambda
            [un_down,ind_down] = max(working_Un(:,1:N),[],2);
            UnDiff_down = un_down - working_Un(:,N);

            pastThreshold_down = UnDiff_down > lambda; % channels that are pass the threshold
            working_Un(increasing_t < ind_down & pastThreshold_down) = NaN; % delete points before inflection so upwards threshold can be checked

            % set imprint to be true between inflection point and now for every channel past theshold upwards
            imprint(increasing_t >= ind & increasing_t <= N & pastThreshold_up) = true; 
            % set imprint to be false between inflection point and now for every channel past theshold downwards
            imprint(increasing_t >= ind_down & increasing_t <= N & pastThreshold_down) = false;
        end

        nchr=sum(sum(imprint,2)>=1);%get number of channels ever in imprint
        onset_time = find(sum(imprint),1, 'first');
        if isempty(onset_time)
            onset = [];
            onset_time = NaN;
        else
            onset = sum(imprint(:,onset_time+0:7),2)>=1;
        end

        % write output
        cell_imprint(s,:) = table(data_tbl.segment_id(s), {imprint});
        cell_madscores(s,:) = table({mahal_mad_mat},{Un});
        cell_pre_mahal_mat{s} = pre_mahal_mat;
        cell_pre_features_mad{s} = pre_mahal_mad_mat;
        cell_t{s}=tw(ictal_ids);

        tbl_imprint_out(s,3:5) = table((onset_time-1)/8, {onset}, nchr);
    end
end