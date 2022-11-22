% Loads data from experiments
%
% Settings:
% (what to load)
% load_parts.(cam/ephys)
%
% (ephys)
% site (if multiple probes)
% recording (if multiple on/off recordings within probe)

%% Initialize paths
myPaths;

%% Display progress or not
if ~exist('verbose', 'var')
    verbose = false;
end

%% Debug mode: plot a buc=nch of useful stuff
if ~exist('debug', 'var')
    debug = false;
end

%% Define what to load

% Site (multiple probes) is optional
if ~exist('site', 'var')
    site = [];
end

% Recording (multiple recordings on one probe) is optional
if ~exist('recording', 'var')
    recording = [];
end

% If nothing specified, load everything (but not LFP)
if ~exist('load_parts', 'var')
    load_parts.cam = true;
    load_parts.ephys = true;
else
    % If only some things specified, don't load others
    if ~isfield(load_parts, 'cam')
        load_parts.cam = false;
    end
    if ~isfield(load_parts, 'ephys')
        load_parts.ephys = false;
    end
end

%% Load timeline and associated inputs

[timeline_filename, timeline_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'timeline');
if ~timeline_exists
    error([animal, ' ', day, ' ', num2str(experiment), ': no timeline']);
end

if timeline_exists
    if verbose
        disp('Loading timeline...');
    end

    load(timeline_filename);

    % Get camera times
    cam_name = 'pcoExposure';
    timeline_cam_idx = strcmp({Timeline.hw.inputs.name}, cam_name);

    cam_expose_starts = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1, timeline_cam_idx) <= 2 & ...
        Timeline.rawDAQData(2:end, timeline_cam_idx) > 2)+1);
    cam_expose_stops = Timeline.rawDAQTimestamps( ...
        find(Timeline.rawDAQData(1:end-1, timeline_cam_idx) >= 2 & ...
        Timeline.rawDAQData(2:end, timeline_cam_idx) < 2)+1);

    cam_time = cam_expose_starts;
    cam_expose_times = cam_expose_stops - cam_expose_starts;

    % Get acqLive signal
    acqLive_name = 'acqLive';
    acqLive_idx = strcmp({Timeline.hw.inputs.name}, acqLive_name);
    thresh = max(Timeline.rawDAQData(:, acqLive_idx)) / 2;
    acqLive_trace = Timeline.rawDAQData(:, acqLive_idx) > thresh;
    acqLive_timeline = Timeline.rawDAQTimestamps( ...
        [find(acqLive_trace, 1), find(acqLive_trace, 1, 'last') + 1]);

    % Get wheel position
    rotaryEncoder_idx = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
    % (this is a very strange hack to overcome a problem in the rotary
    % encoder that's known in the lab and was put on the wiki)
    wheel_position = Timeline.rawDAQData(:, rotaryEncoder_idx);
    wheel_position(wheel_position > 2^31) = wheel_position(wheel_position > 2^31) - 2^32;

    % Get wheel velocity by smoothing the wheel trace and taking deriv
    wheel_smooth_t = 0.05; % seconds
    wheel_smooth_samples = wheel_smooth_t / Timeline.hw.samplingInterval;
    wheel_velocity = interp1(conv(Timeline.rawDAQTimestamps, [1, 1]/2, 'valid'), ...
        diff(smooth(wheel_position, wheel_smooth_samples)), Timeline.rawDAQTimestamps)';

    % Get whether stim was flickering - this is only for when widefiled is
    % combined
    stimScreen_idx = strcmp({Timeline.hw.inputs.name}, 'stimScreen');
    if any(stimScreen_idx)
        stimScreen_flicker = max(Timeline.rawDAQData(:, stimScreen_idx)) - ...
            min(Timeline.rawDAQData(:, stimScreen_idx)) > 2;
    end

    % Get photodiode flips (compensate for screen flicker)
    photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
    stimScreen_on = Timeline.rawDAQData(:, photodiode_idx) > 0.2;
    stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);
    photodiode_thresh = 2; % old: max(Timeline.rawDAQData(:,photodiode_idx))/2
    photodiode_trace = Timeline.rawDAQData(stimScreen_on, photodiode_idx) > photodiode_thresh;

    photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, ...
        photodiode_idx), 3);
    photodiode_diff_thresh = range(Timeline.rawDAQData(:, photodiode_idx)) * 0.2;
    photodiode_diff_t = 20; % time (in ms) to get delayed differential
    photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
    photodiode_diff_filt = [1, zeros(1, photodiode_diff_samples), -1];
    photodiode_diff_conv = abs(conv(photodiode_trace_medfilt, photodiode_diff_filt, 'valid'));
    photodiode_trace_diff = photodiode_diff_conv > ...
        photodiode_diff_thresh;
    photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
        photodiode_trace_diff(2:end)) + photodiode_diff_samples + 1;
    photodiode_flip_times = stimScreen_on_t(photodiode_flip)';


    % Get flipper signal (this was added late, might not be present)
    flipper_name = 'flipper';
    flipper_idx = strcmp({Timeline.hw.inputs.name}, flipper_name);
    flipper_thresh = 2; % TTL threshold
    flipper_trace = Timeline.rawDAQData(:, flipper_idx) > flipper_thresh;
    flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
        (flipper_trace(1:end-1) & ~flipper_trace(2:end))) + 1;
    flipper_flip_times_timeline = Timeline.rawDAQTimestamps(flipper_flip)';

    if debug
        %size(Timeline.rawDAQData)
        samples_to_plot = 29000;
        figure('Color', 'white');
        clf
        subplot(311)
        title('photodiode')
        hold on;
        scatter(photodiode_flip(find(photodiode_flip <= samples_to_plot)), ones(size(find(photodiode_flip <= samples_to_plot), 1), 1), 8, 'filled')
        plot(Timeline.rawDAQData(1:samples_to_plot, photodiode_idx))
        plot(photodiode_diff_conv(1:samples_to_plot))
        xlabel('time (in samples)')


        subplot(312)
        title('timeline flipper')
        hold on;
        plot(Timeline.rawDAQData(1:samples_to_plot, flipper_idx))
        scatter(flipper_flip(find(flipper_flip <= samples_to_plot)), ones(size(find(flipper_flip <= samples_to_plot), 1), 1), 8, 'filled')
        plot(flipper_trace(1:samples_to_plot))
        xlabel('time (in samples)')

        subplot(313)   
        title('acq live')
        hold on;
        plot(Timeline.rawDAQData(1:samples_to_plot, acqLive_idx))
        xlabel('time (in samples)')
    end

end

%% Load task/behavior

% Load the block
[block_filename, block_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'block');

if block_exists

    if verbose
        disp('Loading block file...');
    end

    load(block_filename);

    signals_events = block.events;

    % For task behavior data, if reward information exists, use that to align signals/timeline
    % (bad now because manual reward possible - use flipper in future)
    if exist('Timeline', 'var') && isfield(block.outputs, 'rewardTimes')
        reward_t_block = block.outputs.rewardTimes(block.outputs.rewardValues > 0);

        timeline_reward_idx = strcmp({Timeline.hw.inputs.name}, 'rewardEcho');
        reward_thresh = max(Timeline.rawDAQData(:, timeline_reward_idx)) / 2;
        reward_trace = Timeline.rawDAQData(:, timeline_reward_idx) > reward_thresh;
        reward_t_timeline = Timeline.rawDAQTimestamps(find(reward_trace(2:end) & ~reward_trace(1:end-1))+1);

        if length(reward_t_block) == length(reward_t_timeline) && ...
                ~isempty(reward_t_block)
            % If same number of timeline/block rewards, use those for sync
            timeline2block = reward_t_timeline;
            block2timeline = reward_t_block;

        elseif length(reward_t_block) ~= length(reward_t_timeline) && ...
                ~isempty(reward_t_block)
            % If there's a different number of block and timeline rewards (aka
            % manual rewards were given), try to fix this
            % (this is really inelegant but I think works - find the most
            % common offset between block/timeline rewards)
            reward_t_offset = bsxfun(@minus, reward_t_block', reward_t_timeline);
            blunt_reward_offset = mode(round(reward_t_offset(:)*10)) / 10;
            reward_t_offset_shift = reward_t_offset - blunt_reward_offset;
            t_offset_tolerance = 0.1;
            reward_t_offset_binary = abs(reward_t_offset_shift) < t_offset_tolerance;
            if all(sum(reward_t_offset_binary, 2) == 1)
                % one timeline reward for each block reward, you're good
                % (eliminate the timeline rewards with no match)
                manual_timeline_rewards = sum(reward_t_offset_binary, 1) == 0;
                reward_t_timeline(manual_timeline_rewards) = [];
                warning('Manual rewards included - removed successfully');

                timeline2block = reward_t_timeline;
                block2timeline = reward_t_block;
            else
                % otherwise, you're in trouble
                error('Manual rewards included - couldn''t match to block');
            end
        
    elseif isempty(reward_t_block)
        % If no rewards: probably much less robust, but use stim on
        % times and photodiode flips
        warning('No rewards, aligning block with estimated stimOn times');
        stim_t_offset = block.events.stimOnTimes' - photodiode_flip_times';
        blunt_stim_offset = mode(round(stim_t_offset(:)*10)) / 10;
        stim_t_offset_shift = stim_t_offset - blunt_stim_offset;
        t_offset_tolerance = 0.05;
        stim_t_offset_binary = abs(stim_t_offset_shift) < t_offset_tolerance;
        stim_tl_block_matched = sum(stim_t_offset_binary, 2) == 1;
        if nanmean(stim_tl_block_matched) > 0.9
            % (if 90% stim matched, use this)
            [~, estimated_stimOn_idx] = max(stim_t_offset_binary, [], 2);

            timeline2block = photodiode_flip_times(estimated_stimOn_idx(stim_tl_block_matched));
            block2timeline = block.events.stimOnTimes(stim_tl_block_matched);
        else
            % (if >10% unmatched, don't use)
            error('Attempted stim on alignment - bad match');
        end


        % Go through all block events and convert to timeline time
        % (uses reward as reference)
        block_fieldnames = fieldnames(block.events);
        block_values_idx = cellfun(@(x) ~isempty(x), strfind(block_fieldnames, 'Values'));
        block_times_idx = cellfun(@(x) ~isempty(x), strfind(block_fieldnames, 'Times'));
        for curr_times = find(block_times_idx)'
            if isempty(signals_events.(block_fieldnames{curr_times}))
                % skip if empty
                continue
            end
            signals_events.(block_fieldnames{curr_times}) = ...
                interp1(block2timeline, timeline2block, block.events.(block_fieldnames{curr_times}), 'linear', 'extrap');
        end
        end
    end

    % SPECIFIC TO PROTOCOL
    JF_pull_signals;
    if debug
    end
end

%% Load face/eyecam and processing

% Don't load if no timeline
if exist('Timeline', 'var') && load_parts.cam

    % Get cam sync from timeline
    camSync_idx = strcmp({Timeline.hw.inputs.name}, 'camSync');
    camSync_thresh = max(Timeline.rawDAQData(:, camSync_idx)) / 2;
    camSync = Timeline.rawDAQData(:, camSync_idx) > camSync_thresh;
    camSync_flip = find((camSync(1:end-1) ~= camSync(2:end))) + 1;
    if length(camSync_flip) ~= 4
        error('camSync flip number ~= 4')
    end

    % EYECAM
    [eyecam_dir, eyecam_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'eyecam');

    if eyecam_exists
        if verbose
            disp('Loading eyecam...');
        end

        % Load camera processed data
        [eyecam_processed_filename, eyecam_processed_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'eyecam_processed');
        if eyecam_processed_exists
            eyecam = load(eyecam_processed_filename);
        end

        % Get camera times
        eyecam_fn = AP_cortexlab_filenameJF(animal, day, experiment, 'eyecam');
        eyecam_dir = fileparts(eyecam_fn);
        eyecam_t_savefile = [eyecam_dir, filesep, 'eyecam_t.mat'];

        if exist(eyecam_fn, 'file') && ~exist(eyecam_t_savefile, 'file')
            % Get facecam strobes
            eyeCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'eyeCameraStrobe') | ...
                strcmp({Timeline.hw.inputs.name}, 'eyeCamStrobe');
            eyeCamStrobe_thresh = max(Timeline.rawDAQData(:, eyeCamStrobe_idx)) / 5;
            eyeCamStrobe = Timeline.rawDAQData(:, eyeCamStrobe_idx) > eyeCamStrobe_thresh;
            eyeCamStrobe_up = find((~eyeCamStrobe(1:end-1) & eyeCamStrobe(2:end))) + 1;
            eyeCamStrobe_up_t = Timeline.rawDAQTimestamps(eyeCamStrobe_up);

            % Get sync times for cameras (or load if already done)
            [eyecam_sync_frames, n_eyecam_frames] = AP_get_cam_sync_framesJF(eyecam_fn);

            if ~isempty(eyecam_sync_frames)
                % Get the closest cam strobe to sync start, find offset and frame idx
                [~, eyecam_strobe_sync] = min(abs(camSync_flip(1)-eyeCamStrobe_up));
                eyecam_frame_offset = eyecam_sync_frames(1) - eyecam_strobe_sync;
                eyecam_frame_idx = [1:length(eyeCamStrobe_up)] + eyecam_frame_offset;

                % Check that the number of frames between synchs matches
                % video and timeline
                n_eyecam_frames_syncd_movie = diff(eyecam_sync_frames) + 1;
                [~, eyecam_strobe_sync_end] = min(abs(camSync_flip(3)-eyeCamStrobe_up));
                n_eyecam_frames_syncd_timeline = eyecam_strobe_sync_end - eyecam_strobe_sync;
                if abs(n_eyecam_frames_syncd_movie-n_eyecam_frames_syncd_timeline) > 2
                    warning('Eyecam: different n frames video vs timeline');
                end

                % Get times of cam frames in timeline
                eyecam_t = nan(n_eyecam_frames, 1);
                eyecam_t(eyecam_frame_idx(eyecam_frame_idx > 0)) = eyeCamStrobe_up_t(eyecam_frame_idx > 0);

                save(eyecam_t_savefile, 'eyecam_t');
            end
        elseif exist(eyecam_fn, 'file') && exist(eyecam_t_savefile, 'file')
            load(eyecam_t_savefile);
        end

    end

    % FACECAM
    [facecam_dir, facecam_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'facecam');

    if facecam_exists
        if verbose;
            disp('Loading facecam...');
        end

        % Get camera times
        facecam_fn = AP_cortexlab_filenameJF(animal, day, experiment, 'facecam');
        facecam_dir = fileparts(facecam_fn);
        facecam_t_savefile = [facecam_dir, filesep, 'facecam_t.mat'];

        if exist(facecam_fn, 'file') && ~exist(facecam_t_savefile, 'file')
            % Get facecam strobes
            faceCamStrobe_idx = strcmp({Timeline.hw.inputs.name}, 'faceCamStrobe');
            faceCamStrobe_thresh = max(Timeline.rawDAQData(:, faceCamStrobe_idx)) / 5;
            faceCamStrobe = Timeline.rawDAQData(:, faceCamStrobe_idx) > faceCamStrobe_thresh;
            faceCamStrobe_up = find((~faceCamStrobe(1:end-1) & faceCamStrobe(2:end))) + 1;
            faceCamStrobe_up_t = Timeline.rawDAQTimestamps(faceCamStrobe_up);

            % Get sync times for cameras (or load if already done)
            [facecam_sync_frames, n_facecam_frames] = AP_get_cam_sync_framesJF(facecam_fn);

            if ~isempty(facecam_sync_frames)
                % Get the closest cam strobe to sync start, find offset and frame idx
                [~, facecam_strobe_sync] = min(abs(camSync_flip(1)-faceCamStrobe_up));
                facecam_frame_offset = facecam_sync_frames(1) - facecam_strobe_sync;
                facecam_frame_idx = [1:length(faceCamStrobe_up)] + facecam_frame_offset;

                % Check that the number of frames between syncs matches
                % video and timeline
                n_facecam_frames_syncd_movie = diff(facecam_sync_frames) + 1;
                [~, facecam_strobe_sync_end] = min(abs(camSync_flip(3)-faceCamStrobe_up));
                n_facecam_frames_syncd_timeline = facecam_strobe_sync_end - facecam_strobe_sync;
                if abs(n_facecam_frames_syncd_movie-n_facecam_frames_syncd_timeline) > 2
                    warning('Facecam: different n frames video vs timeline');
                end

                % Get times of cam frames in timeline
                facecam_t = nan(n_facecam_frames, 1);
                facecam_t(facecam_frame_idx(facecam_frame_idx > 0)) = faceCamStrobe_up_t(facecam_frame_idx > 0);

                save(facecam_t_savefile, 'facecam_t');
            end
        elseif exist(facecam_fn, 'file') && exist(facecam_t_savefile, 'file')
            load(facecam_t_savefile);
        end

        % (old/unused: etGUI and facemap)
        [facecam_processed_filename, facecam_processed_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'facecam_processed');
        if facecam_processed_exists
            facecam = load(facecam_processed_filename);
        end

        % (output from AP_mouse_movie_movement)
        [facecam_movement_filename, facecam_movement_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'facecam_movement');
        if facecam_movement_exists
            load(facecam_movement_filename);
        end

    end

end

%% Load ephys data (single long recording)

% Pick kilosort version (2 by default, 1 old if selected)
if ~exist('kilosort_version', 'var') || kilosort_version == 2
    [ephys_path, ephys_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys', site, recording);
elseif exist('kilosort_version', 'var') && kilosort_version == 1
    [ephys_path, ephys_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_ks1', site, recording);
end

if ephys_exists && load_parts.ephys

    if verbose;
        disp('Loading ephys...');
    end

    % These are the digital channels going into the FPGA
    photodiode_sync_idx = 1;
    acqLive_sync_idx = 2;
    led_sync_idx = 3;
    flipper_sync_idx = 4;

    % Load phy sorting if it exists
    % (old = cluster_groups.csv, new = cluster_group.tsv because fuck me)
    cluster_filepattern = [ephys_path, filesep, 'cluster_group*'];
    cluster_filedir = dir(cluster_filepattern);
    if ~isempty(cluster_filedir)
        cluster_filename = [ephys_path, filesep, cluster_filedir.name];
        fid = fopen(cluster_filename);
        cluster_groups = textscan(fid, '%d%s', 'HeaderLines', 1);
        cg = cluster_groups;
        fclose(fid);
    end

    % Load sync/photodiode
    if isempty(dir([ephys_path, filesep, 'sync.mat']))
        error('No SYNC - extract sync before continuing')
    end
    load(([ephys_path, filesep, 'sync.mat']));

    if debug
    end

    % Read header information
    if ~isSpikeGlx
        header_path = [ephys_path, filesep, 'dat_params.txt'];
        header_fid = fopen(header_path);
        header_info = textscan(header_fid, '%s %s', 'delimiter', {' = '});
        fclose(header_fid);

        header = struct;
        for i = 1:length(header_info{1})
            header.(header_info{1}{i}) = header_info{2}{i};
        end

        % Load spike data
        if isfield(header, 'sample_rate')
            ephys_sample_rate = str2num(header.sample_rate);
        elseif isfield(header, 'ap_sample_rate')
            ephys_sample_rate = str2num(header.ap_sample_rate);
        end
    else
        header.n_channels = 385;
        ephys_sample_rate = 30000;
        header.lfp_sample_rate = 30000;
        header.filter_cutoff = 300; %Hz, defined in JF_computeLFP

    end
    spike_times = double(readNPY([ephys_path, filesep, 'spike_times.npy'])) ./ ephys_sample_rate;
    spike_templates_0idx = readNPY([ephys_path, filesep, 'spike_templates.npy']);
    templates_whitened = readNPY([ephys_path, filesep, 'templates.npy']);
    channel_positions = readNPY([ephys_path, filesep, 'channel_positions.npy']);
    channel_map = readNPY([ephys_path, filesep, 'channel_map.npy']);
    winv = readNPY([ephys_path, filesep, 'whitening_mat_inv.npy']);
    template_amplitudes = readNPY([ephys_path, filesep, 'amplitudes.npy']);

    % Default channel map/positions are from end: make from surface
    % (hardcode this: kilosort2 drops channels)
    if isSpikeGlx
        max_depth = 2880;
        if any(max_depth-channel_positions(:, 2) < 0) %1.0
            max_depth = 3840;
            channel_positions(:, 2) = max_depth - channel_positions(:, 2);
        else
            max_depth = 2880;
            channel_positions(:, 2) = max_depth - channel_positions(:, 2);
        end

    else
        max_depth = 3840;
        channel_positions(:, 2) = max_depth - channel_positions(:, 2); % 0 = tip in NP1s, 3840 = top, reorder here
    end


    % Unwhiten templates
    templates = zeros(size(templates_whitened));
    for t = 1:size(templates_whitened, 1)
        templates(t, :, :) = squeeze(templates_whitened(t, :, :)) * winv;
    end

    % Get the waveform of all templates (channel with largest amplitude)
    [~, max_site] = max(max(abs(templates), [], 2), [], 3);
    templates_max = nan(size(templates, 1), size(templates, 2));
    for curr_template = 1:size(templates, 1)
        templates_max(curr_template, :) = ...
            templates(curr_template, :, max_site(curr_template));
    end
    waveforms = templates_max;

    % Get depth of each template
    % (get min-max range for each channel)
    template_chan_amp = squeeze(range(templates, 2));
    % (zero-out low amplitude channels)
    template_chan_amp_thresh = max(template_chan_amp, [], 2) * 0.5;
    template_chan_amp_overthresh = template_chan_amp .* (template_chan_amp >= template_chan_amp_thresh);
    % (get center-of-mass on thresholded channel amplitudes)
    template_depths = sum(template_chan_amp_overthresh.*channel_positions(:, 2)', 2) ./ sum(template_chan_amp_overthresh, 2);
    template_xdepths = channel_positions(max_site, 1);

    % Get the depth of each spike (templates are zero-indexed)
    spike_depths = template_depths(spike_templates_0idx+1);
    spike_xdepths = template_xdepths(spike_templates_0idx+1);

    % Get trough-to-peak time for each template
    templates_max_signfix = bsxfun(@times, templates_max, ...
        sign(abs(min(templates_max, [], 2))-abs(max(templates_max, [], 2))));

    [~, waveform_trough] = min(templates_max, [], 2);
    [~, waveform_peak_rel] = arrayfun(@(x) ...
        max(templates_max(x, waveform_trough(x):end), [], 2), ...
        transpose(1:size(templates_max, 1)));
    waveform_peak = waveform_peak_rel + waveform_trough;

    templateDuration = waveform_peak - waveform_trough;
    templateDuration_us = (templateDuration / ephys_sample_rate) * 1e6;

    % Get sync points for alignment

    % Get experiment index by finding numbered folders
    protocols_list = AP_list_experimentsJF(animal, day);
    experiment_idx = experiment == [protocols_list.experiment];

    % load sync and align 
    [spike_times_timeline, bad_flipper] = JF_align_ephys_to_timeline(animal, day, isSpikeGlx, flipper_flip_times_timeline, ...
    ephys_path, flipper_sync_idx, experiment_idx, acqLive_sync_idx, spike_times);
 

        %         co = robustfit(sync_ephys, sync_timeline);
        %          spike_times_timeline = spike_times * co(2) + co(1);
    
    % Get "good" templates from labels
    if exist('cluster_groups', 'var') && loadClusters
        % If there's a manual classification


        % Check that all used spike templates have a label
        spike_templates_0idx_unique = unique(spike_templates_0idx);
        if ~all(ismember(spike_templates_0idx_unique, uint32(cluster_groups{1}))) || ...
                ~all(ismember(cluster_groups{2}, {'good', 'mua', 'noise'}))
            disp('Removing manually labelled noise units...');
            good_templates_idx = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'unsorted') | strcmp(cluster_groups{2}, 'mua') | ...
                strcmp(cluster_groups{2}, 'good')));
            template_label_mua = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'mua')));
            template_label_good = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'good')));
            template_labelM = good_templates_idx(ismember(good_templates_idx, template_label_mua));
            template_labelG = good_templates_idx(ismember(good_templates_idx, template_label_good));

            [good_templates, ~] = ismember(0:size(templates, 1)-1, good_templates_idx);
            [good_templatesG, ~] = ismember(0:size(templates, 1)-1, template_labelG);
            [good_templatesM, ~] = ismember(0:size(templates, 1)-1, template_labelM);


        else
            disp('Keeping manually labelled good units...');
            good_templates_idx = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'good') | strcmp(cluster_groups{2}, 'mua')));
            template_label_mua = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'mua')));
            template_label_good = uint32(cluster_groups{1}( ...
                strcmp(cluster_groups{2}, 'good')));
            template_labelM = good_templates_idx(ismember(good_templates_idx, template_label_mua));
            template_labelG = good_templates_idx(ismember(good_templates_idx, template_label_good));

            [good_templates, ~] = ismember(0:size(templates, 1)-1, good_templates_idx);
            [good_templatesG, ~] = ismember(0:size(templates, 1)-1, template_labelG);
            [good_templatesM, ~] = ismember(0:size(templates, 1)-1, template_labelM);

            template_label = good_templatesG .* 1 + (good_templatesM) .* 2;
        end


    elseif exist('unitType', 'var')
        % If no manual but qualityMetrics are available
        if verbose;
            disp('Keeping quality metrics good units...');
        end

        % Load triage labels

        %triage_good_templates = goodUnits;

        good_templates = ...
            unitType == 1;
        good_templates_idx = find(unitType == 1) - 1;

        if exist('locationKeep', 'var')
            if verbose;
                disp('Keeping location data...');
            end
            myPaths;
            allenAt = loadStructureTreeJF([allenAtlasPath, filesep, 'allenCCF/structure_tree_safe_2017.csv']);
            probeccf = AP_cortexlab_filenameJF(animal, [], [], 'histo', [], []);
            load(probeccf)
            this_ccf = probe_ccf(probes(iDataset));
            theseLocations = allenAt.acronym(this_ccf.trajectory_areas);
            theseLocationsInterest = contains(theseLocations, locationKeep);
            theseDepths = this_ccf.probe_depths(theseLocationsInterest);
            template_exists = ismember(1:max(spikeTemplates), unique(spikeTemplates));
            theseTemplates = template_depths(template_exists) >= min(theseDepths) & template_depths(template_exists) <= max(theseDepths); %correct depth units
            good_templates = goodUnits & theseTemplates';
            good_templates_idx = find(good_templates) - 1;
        end
    else
        % If no cluster groups at all, keep all
        warning([animal, ' ', day, ' - no cluster groups']);
        if verbose
            disp('No manual labeling, keeping all and re-indexing');
        end
        good_templates_idx = unique(spike_templates_0idx);
        good_templates = ismember(0:size(templates, 1)-1, good_templates_idx);
    end
    % Throw out all non-good template data
    templates = templates(good_templates, :, :);
    template_depths = template_depths(good_templates);
    template_xdepths = template_xdepths(good_templates);
    waveforms = waveforms(good_templates, :);
    templateDuration = templateDuration(good_templates);
    templateDuration_us = templateDuration_us(good_templates);
    %template_label = template_label(good_templates);
    % Throw out all non-good spike data
    good_spike_idx = ismember(spike_templates_0idx, good_templates_idx);
    spike_times = spike_times(good_spike_idx);
    spike_times_full = spike_times_timeline;
    spike_templates_full = spike_templates_0idx + 1;
    spike_templates_0idx = spike_templates_0idx(good_spike_idx);
    template_amplitudes = template_amplitudes(good_spike_idx);
    spike_depths = spike_depths(good_spike_idx);
    spike_xdepths = spike_xdepths(good_spike_idx);
    spike_times_timeline = spike_times_timeline(good_spike_idx);

    % Rename the spike templates according to the remaining templates
    % (and make 1-indexed from 0-indexed)
    new_spike_idx = nan(max(spike_templates_0idx)+1, 1);
    new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
    spike_templates = new_spike_idx(spike_templates_0idx+1);

end

%% Load LFP
% (either single channel full or all channel snippet)

if ephys_exists && load_parts.ephys && exist('lfp_channel', 'var')

    % Get LFP file info

    if ~isSpikeGlx
        n_channels = str2num(header.n_channels);
        [data_path, data_path_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'ephys_dir', site);


        % (get all recordings within site - assume concat at this point)
        lfp_recordings = dir([data_path, filesep, 'experiment*']);
        lfp_filenames = cellfun(@(x) ...
            [data_path, filesep, x, filesep, 'recording1/continuous/Neuropix-3a-100.1/continuous.dat'], ...
            {lfp_recordings.name}, 'uni', false);

        % Get LFP properties
        % (NOTE: LFP channel map is different from kilosort channel map because
        % kilosort2 drops channels without spikes)
        channel_map_fn = [dropboxPath, filesep, 'MATLAB/onPaths/JF_Scripts_CortexLab/kilosort/forPRBimecP3opt3.mat'];
        channel_map_full = load(channel_map_fn);
        max_depth = 3840;
        lfp_channel_positions = max_depth - channel_map_full.ycoords;
        lfp_channel_xpositions = channel_map_full.xcoords;
        lfp_sample_rate = str2num(header.lfp_sample_rate);
        lfp_cutoff = str2num(header.filter_cutoff);

        % Memory map LFP
        n_bytes = 2; % LFP = int16 = 2 bytes
        n_lfp_samples = nan(size(lfp_filenames));
        lfp_memmap = cell(size(lfp_filenames));
        for curr_lfp_filename = 1:length(lfp_filenames)
            lfp_fileinfo = dir(lfp_filenames{curr_lfp_filename});
            n_lfp_samples(curr_lfp_filename) = lfp_fileinfo.bytes / n_bytes / n_channels;
            lfp_memmap{curr_lfp_filename} = ...
                memmapfile(lfp_filenames{curr_lfp_filename}, ...
                'Format', {'int16', [n_channels, n_lfp_samples(curr_lfp_filename)], 'lfp'});
        end
    else
        [filename, file_exists] = AP_cortexlab_filenameJF(animal, day, experiment, 'mainfolder', site, recording);
        mainFolder = filename{1}(1:end - 1);
        lfpDir = dir([mainFolder, filesep, 'lfp', filesep, 'lfp.mat']);
        if isempty(lfpDir) && loadLFP %only compute LFP if not already saved on disk
            n_channels = header.n_channels;
            lfp = JF_get_NPX2_LFP(ephys_path);
            mkdir([mainFolder, filesep, 'lfp'])
            save([mainFolder, filesep, 'lfp', filesep, 'lfp.mat'], 'lfp', '-v7.3') %save lfp
        elseif loadLFP
            load([mainFolder, filesep, 'lfp', filesep, 'lfp.mat'])
        else
            lfp = NaN;
            lfp_sample_rate = 30000 / 1000;
            binFile = dir([ephys_path(1:end-15), ephys_path(end-4:end), filesep, '*ap.meta']);
            if isempty(binFile)
                binFile = dir([ephys_path(1:end-15), ephys_path(end-4:end), filesep, 'experiment1/*ap.meta']);
            end

            meta = ReadMeta_GLX(binFile.name, binFile.folder);
            ephysAPfile = [ephys_path(1:end-15), ephys_path(end-4:end)];
            %max_depth = 0;%already ordered
            if contains(meta.imRoFile, 'NPtype24_hStripe_botRow0_ref0.imro') %2.0 4 shank, bottom stripe
                chanMapData = load([dropboxPath, filesep, 'MATLAB/onPaths/JF_Scripts_CortexLab/kilosort/chanMapNP2_4Shank_bottRow_flipper.mat']);
                lfp_channel_positions = chanMapData.ycoords;
                lfp_channel_xpositions = chanMapData.xcoords;
            else %1.0 bottom row
                chanMapData = load([dropboxPath, filesep, 'MATLAB/onPaths/JF_Scripts_CortexLab/kilosort/chanMapNP2_1Shank_flipper.mat']);
                lfp_channel_positions = chanMapData.ycoords;
                lfp_channel_xpositions = chanMapData.xcoords;
            end

        end
    end
    if isnumeric(lfp_channel) && loadLFP

        % Load LFP of whole current experiment from one channel
        if verbose;
            disp(['Loading LFP (channel ', num2str(lfp_channel), ')...']);
        end;

        % Load single LFP channel within experiment bounds
        % (treat as concatenated if multiple files)
        lfp_load_start = round((lfp_sample_rate * sync_ephys(1)));
        lfp_load_stop = round((lfp_sample_rate * sync_ephys(end)));

        lfp_concat_length = [0, cumsum(n_lfp_samples)];
        lfp_load_file = find(lfp_load_start < lfp_concat_length, 1) - 1;
        lfp_load_start_rel = lfp_load_start - lfp_concat_length(lfp_load_file);
        lfp_load_stop_rel = lfp_load_stop - lfp_concat_length(lfp_load_file);

        lfp = double(lfp_memmap{lfp_load_file}.Data.lfp(lfp_channel, lfp_load_start_rel:lfp_load_stop_rel));

        % Get LFP times and convert to timeline time
        lfp_load_start_t = lfp_load_start / lfp_sample_rate;
        lfp_t = [0:size(lfp, 2) - 1] / lfp_sample_rate + lfp_load_start_t;
        lfp_t_timeline = interp1(sync_ephys, sync_timeline, lfp_t, 'linear', 'extrap');

        %%% Remove light artifact
        if verbose;
            disp('Cleaning LFP...');
        end;

        % Get light times (assume blue/violet alternate)
        light_t_timeline = interp1(sync_ephys, sync_timeline, sync(led_sync_idx).timestamps, 'linear', 'extrap');
        use_light_times = light_t_timeline >= lfp_t_timeline(1) & light_t_timeline <= lfp_t_timeline(end);
        light_on = light_t_timeline(sync(led_sync_idx).values == 1 & use_light_times);
        light_off = light_t_timeline(sync(led_sync_idx).values == 0 & use_light_times);

        % (cut uncoupled off/on from start/end)
        if light_off(1) < light_on(1)
            light_off(1) = [];
        end
        if light_on(end) > light_off(end)
            light_on(end) = [];
        end

        blue_on = light_on(1:2:end);
        blue_off = light_off(1:2:end);
        violet_on = light_on(2:2:end);
        violet_off = light_off(2:2:end);

        light_on_mean = mean(light_off-light_on);
        light_off_mean = mean(light_on(2:end)-light_off(1:end-1));
        light_surround_t = [-(light_off_mean / 2):1 / lfp_sample_rate:(light_on_mean + (light_off_mean / 2))];

        % Pull out LFP around light on
        use_blue_on = blue_on >= lfp_t_timeline(1) & blue_on <= lfp_t_timeline(end);
        blue_on_pull_t = blue_on(use_blue_on) + light_surround_t;
        blue_on_lfp = interp1(lfp_t_timeline, lfp', blue_on_pull_t);

        use_violet_on = violet_on >= lfp_t_timeline(1) & violet_on <= lfp_t_timeline(end);
        violet_on_pull_t = violet_on(use_violet_on) + light_surround_t;
        violet_on_lfp = interp1(lfp_t_timeline, lfp', violet_on_pull_t);

        % Subtract baseline
        baseline_t = find(light_surround_t < 0, 1, 'last');
        blue_on_lfp_baselinesub = blue_on_lfp - blue_on_lfp(:, baseline_t, :);
        violet_on_lfp_baselinesub = violet_on_lfp - violet_on_lfp(:, baseline_t, :);

        % Get rolling median (allow light artifact to change slightly)
        n_light = 500;
        blue_on_lfp_baselinesub_med = movmedian(blue_on_lfp_baselinesub, n_light, 1);
        violet_on_lfp_baselinesub_med = movmedian(violet_on_lfp_baselinesub, n_light, 1);

        % Interpolate out the artifact to remove
        n_lfp_channels = size(lfp, 1);
        blue_light_remove = interp1( ...
            reshape(permute(blue_on_pull_t, [2, 1]), [], 1), ...
            reshape(permute(blue_on_lfp_baselinesub_med, [2, 1, 3]), [], n_lfp_channels), ...
            reshape(lfp_t_timeline, [], 1))';
        violet_light_remove = interp1( ...
            reshape(permute(violet_on_pull_t, [2, 1]), [], 1), ...
            reshape(permute(violet_on_lfp_baselinesub_med, [2, 1, 3]), [], n_lfp_channels), ...
            reshape(lfp_t_timeline, [], 1))';

        % Zero-out any NaNs (e.g. remove nothing)
        blue_light_remove(isnan(blue_light_remove)) = 0;
        violet_light_remove(isnan(violet_light_remove)) = 0;

        % Remove the artifact
        lfp_lightfix = lfp - (blue_light_remove + violet_light_remove);

        % NOT DOING THIS: IS THIS NECESSARY? TOO MUCH MEMORY
        %     % (low-pass filter: sometimes bunch of junk at high freq?)
        %     freqCutoff = 300; % Hz
        %     [b100s, a100s] = butter(2,freqCutoff/(lfp_sample_rate/2),'low');
        %     lfp_lightfix = single(filtfilt(b100s,a100s,double(lfp_lightfix)')');

    elseif strcmp(lfp_channel, 'all') && loadLFP

        % Load short LFP segment (from start = no light) from all channels
        if verbose
            disp('Loading LFP (all channels snippet)...');
        end
        if ~isSpikeGlx
            % Choose snippet of recording time before first experiment (no light)
            t_load = 10; % time to load (in seconds)
            t_load_pre_exp = 1; % time before first experiment to load up to

            experiment_ephys_starts = sync(acqLive_sync_idx).timestamps(sync(acqLive_sync_idx).values == 1);
            lfp_load_start = round((lfp_sample_rate * (experiment_ephys_starts(1) - t_load_pre_exp - t_load)));
            lfp_load_stop = round((lfp_sample_rate * (experiment_ephys_starts(1) - t_load_pre_exp)));

            % Load all LFP channels in snippet (before recording = first file)
            if lfp_load_stop > size(lfp_memmap, 2)
                lfp = lfp_memmap{1}.Data.lfp(:, 1:lfp_load_stop-lfp_load_start);
            else
                lfp = lfp_memmap{1}.Data.lfp(:, lfp_load_start:lfp_load_stop);
            end
            % Sort LFP so it goes from surface to depth
            [~, lfp_sort_idx] = sort(lfp_channel_positions);
            lfp_channel_positions = lfp_channel_positions(lfp_sort_idx);
            lfp = lfp(lfp_sort_idx, :);

            % Get LFP times and convert to timeline time
            lfp_load_start_t = lfp_load_start / lfp_sample_rate;
            lfp_t = [0:size(lfp, 2) - 1] / lfp_sample_rate + lfp_load_start_t;
            window_length = 2; % in seconds
            window_overlap = 1; % in seconds
            window_length_samples = round(window_length/(1 / lfp_sample_rate));
            window_overlap_samples = round(window_overlap/(1 / lfp_sample_rate));
            [lfp_power, lfp_power_freq] = pwelch(zscore(double(lfp), [], 2)', ...
                window_length_samples, window_overlap_samples, [], lfp_sample_rate);
        else
            lfp_t = [0:size(lfp, 1) - 1] / (30000) * 100;
            lfp_sample_rate = 30000 / 1000;
            binFile = dir([ephys_path(1:end-15), ephys_path(end-4:end), filesep, '*.meta']);
            if isempty(binFile)
                binFile = dir([ephys_path(1:end-15), ephys_path(end-4:end), filesep, 'experiment1/*.meta']);
            end

            meta = ReadMeta_GLX(binFile.name, binFile.folder);
            ephysAPfile = [ephys_path(1:end-15), ephys_path(end-4:end)];
            %max_depth = 0;%already ordered
            if contains(meta.imRoFile, 'NPtype24_hStripe_botRow0_ref0.imro') %2.0 4 shank, bottom stripe
                chanMapData = load([dropboxPath, filesep, 'MATLAB/JF_scripts_CortexLab/kilosort/chanMapNP2_4Shank_bottRow_flipper.mat']);
                lfp_channel_positions = chanMapData.ycoords;
                lfp_channel_xpositions = chanMapData.xcoords;
            else %1.0 bottom row
                chanMapData = load([dropboxPath, filesep, 'MATLAB/JF_scripts_CortexLab/kilosort/chanMapNP2_1Shank_flipper.mat']);
                lfp_channel_positions = chanMapData.ycoords;
                lfp_channel_xpositions = chanMapData.xcoords;
            end
            lfp = permute(lfp, [2, 1]);
            [~, lfp_sort_idx] = sort(lfp_channel_positions);
            lfp_channel_positions = lfp_channel_positions(lfp_sort_idx);
            maxIdx = lfp_sort_idx == max(lfp_sort_idx);
            lfp_sort_idx(maxIdx) = [];
            lfp = lfp(lfp_sort_idx, :);

            window_length = 2; % in seconds
            window_overlap = 1; % in seconds
            window_length_samples = round(window_length/(1 / lfp_sample_rate));
            window_overlap_samples = round(window_overlap/(1 / lfp_sample_rate));
            [lfp_power, lfp_power_freq] = pwelch(zscore(double(lfp), [], 2)', ...
                window_length_samples, window_overlap_samples, [], lfp_sample_rate);

        end


        if verbose
            figure;

            p1 = subplot(1, 2, 1);
            imagesc(lfp_power_freq, lfp_channel_positions(1:end-1), log10(lfp_power'));
            xlabel('Frequency');
            ylabel('Depth (\mum)');
            c = colorbar;
            ylabel(c, 'Log_{10} power');
            xlim([0, 100]);
            colormap(p1, hot);
            title('LFP power');

            p2 = subplot(1, 2, 2);
            imagesc(lfp_channel_positions, lfp_channel_positions, ...
                corrcoef((movmedian(zscore(double(lfp), [], 2), 10, 1) - ...
                nanmedian(zscore(double(lfp), [], 2), 1))'));
            axis image
            colormap(p2, brewermap([], '*RdBu'));
            caxis([-1, 1])
            xlabel('Depth (\mum)');
            ylabel('Depth (\mum)');
            c = colorbar;
            ylabel(c, 'Med. sub. correlation');

        end
    end
end

%% Finished
if verbose
    disp('Finished loading experiment.');
end