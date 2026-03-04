%%% Extract Signal from Axon Attached IC Steps Recordings with EAPs

% Author: Sophie Liebergall
% Updated: 7/16/24
% Purpose: Extract signal from loose seal recordings of axonal APs during
% current step

%% Load in .abf file

%clear all

% Set path to .abf file
% pathFileName = "\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\data\axon_propagation\slice_electrophysiology\PV\ICsteps\7391_1_ICsteps.abf";

% get the path of the .abf file of interest
%     [fileName,path] = uigetfile('*.abf', 'abf files');
%     pathFileName = [path fileName];
%     % make this the default path
%     cd(path)

[data,si,h]=abfload(pathFileName);
% data = data frame configured as <data pts per sweep> by <number of chans> 
%   by <number of sweeps>.
% si = sampling interval in microseconds (us)
%       Note: sampling rate (Hz) = 1e6 / si (us)
% h = information on file (selected header parameters)

% Extract dimensions of the full data matrix
[i, num_chan, num_sweeps]=size(data); 
% i = number of samples in each sweep (i.e. recorded data points)
% num_chan = number of channels
% num_sweeps = number of sweeps in .abf file

% Specify channel of interest
% User may need to change this value depending on channel of interest

soma_v_chan = 1;
soma_i_chan = 2;
axon_i_chan = 3;
axon_v_chan = 4;
    
% Create vector with all indices in file
indices = 1:i;

% Create a vector with all indices converted to time points
time(indices,1) = (indices./(1e6./si));

% Create table to store somatic and axon AP data
AP_data = table();

%% Find Indices of Somatic APs

soma_AP_locs = []; % Create array to store indices of soma APs
current_injection = []; % Create array to store current injection for sweep

% Find where injecting at least 100 pA current
current_indices = find(data(:,soma_i_chan,29) > 100);
sweep_length = 0.6; % in sec

for sweep = 1:size(data,3)

% Get current injection for sweep    
current = data(current_indices(1000), soma_i_chan, sweep) - data(current_indices(1) - 1000, soma_i_chan, sweep);
current_injection = [current_injection; current];

trace = data(:,soma_v_chan,sweep);

% Subtract resting membrane potential for peak finding
trace = trace - mean(trace(1:current_indices(1)-1));

% Set threhold for counting successful soma spike
soma_spike_threshold = 50;
% Set minimum time between peaks (i.e. minimum width of spike)
min_spike_width = 0.002; % 2 ms
% Set min on either side before the signal attains a higher value
soma_spike_prom = 20;

fs = 1e6 / si; % sampling rate
[amps,locs] = findpeaks(trace, fs, 'MinPeakHeight', soma_spike_threshold,...
    'MinPeakDistance', min_spike_width, 'MinPeakProminence', soma_spike_prom);
locs = int64(locs * fs);

% figure
% plot(trace)
% hold on
% %Plot triangles at the peak indices
% plot(x(locs), trace(locs), 'v', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

soma_AP_locs = [soma_AP_locs; {locs}];

end

AP_data.soma_AP_locs = soma_AP_locs; % Add to AP data table
AP_data.current_injection = current_injection; % Add to AP data table

current_injection_start = current_indices(1); % index where I injection start
current_injection_end = current_indices(1)+sweep_length*(1e6 / si); % index where I injection ends

% Sort APs into evoked and ectopic
for row = 1:height(AP_data)

    trace = AP_data.soma_AP_locs{row,1};
    % get APs inside current injection
    evoked = trace(trace > current_injection_start & trace < current_injection_end);
    % get APs outside current injection
    EAPs = trace(trace < current_injection_start | trace > current_injection_end);
    
    AP_data.evoked_soma_AP_locs{row,1} = evoked;
    AP_data.ectopic_soma_AP_locs{row,1} = EAPs;

end

%% Find Indices of Axonal APs

% Create arrays to store data across sweeps
 axon_AP_locs_big = {};
 axon_AP_fidelity_big = {};
 axon_AP_latencies_big = {};
 axon_AP_peaks_big = {};

for sweep = 1:size(data,3) % loop over each sweep in recording

soma_AP_locs = AP_data.soma_AP_locs{sweep,1};
    
num_APs = size(soma_AP_locs,1); % number of APs in trace

% Define region to calculate mean and SD of noise
noise_region = data(1:current_injection_start,axon_i_chan,sweep);

%mean_noise = mean(noise_region); % get mean of noise
sd_noise = std(noise_region); % get SD of noise

% Set threhold for counting successful axon spike
axon_spike_threshold = 3*sd_noise;

% Create arrays to store data (each row is a sweep)
axon_AP_locs = [];
axon_AP_peaks = [];
axon_AP_latencies = [];
axon_AP_fidelity = [];

axon_AP_window = 0.5; % Set window (one sided) to detect axon apike (in ms)
axon_AP_window_i = axon_AP_window/(1e-3*si); % convert to indices

for spike = 1:size(soma_AP_locs,1) % loop over each spike in trace

    soma_AP_peak_i = soma_AP_locs(spike); % get index of soma AP peak
    
    % Select region to search for axon spike
    search_region = data(soma_AP_peak_i-axon_AP_window_i:soma_AP_peak_i+axon_AP_window_i,...
        axon_i_chan,sweep);

    % Get local minimum of search region
    [axon_AP_inward, axon_AP_peak_i] = min(search_region);
    axon_AP_peak_i = axon_AP_peak_i + (soma_AP_peak_i-axon_AP_window_i);
    % Get mean value 1 to 0.5 ms before AP peak
    if axon_AP_peak_i < axon_AP_window_i*2
    baseline = mean(data(1:axon_AP_peak_i-(axon_AP_window_i),...
        axon_i_chan,sweep));
    else
    baseline = mean(data(axon_AP_peak_i-(axon_AP_window_i*2):axon_AP_peak_i-(axon_AP_window_i),...
        axon_i_chan,sweep));
    end
    axon_AP_peak = abs(axon_AP_inward - baseline);

    % Get first local max of outward current after axon spike peak
    width_window = 0.8; % Set width of window to look for max outward current after peak (in ms)
    width_window = width_window/(1e-3*si); % convert to indices
    width_search_region = data(axon_AP_peak_i:axon_AP_peak_i+width_window,axon_i_chan,sweep);
    [axon_AP_outward,outward_max_i] = max(width_search_region);
    
    axon_AP_locs = [axon_AP_locs axon_AP_peak_i];
    axon_AP_peaks = [axon_AP_peaks axon_AP_peak];

    if axon_AP_peak > axon_spike_threshold
        axon_AP_fidelity = [axon_AP_fidelity 1];
    else
        axon_AP_fidelity = [axon_AP_fidelity 0];
    end

    % Calculate latency from somatic spike peak to axon spike peak
    latency = axon_AP_peak_i - soma_AP_peak_i;
    latency = double(latency) * (1e-3*si); % convert to ms
    axon_AP_latencies = [axon_AP_latencies latency];

end
 
 axon_AP_locs_big{sweep, 1} = axon_AP_locs;
 axon_AP_fidelity_big{sweep, 1} = axon_AP_fidelity;
 axon_AP_latencies_big{sweep, 1} = axon_AP_latencies;
 axon_AP_peaks_big{sweep, 1} = axon_AP_peaks;
 

end

AP_data.axon_AP_locs = axon_AP_locs_big;
AP_data.axon_AP_fidelity = axon_AP_fidelity_big;
AP_data.axon_AP_latencies = axon_AP_latencies_big;
AP_data.axon_AP_peaks = axon_AP_peaks_big;

% Sort APs into evoked and ectopic
for row = 1:height(AP_data)

    trace = AP_data.axon_AP_locs{row,1};
    % get APs inside current injection
    evoked_i = (trace > current_injection_start & trace < current_injection_end);
    % get APs outside current injection
    ectopic_i = (trace < current_injection_start | trace > current_injection_end);
    
    evoked_axon_AP_locs = axon_AP_locs_big{row,1};
    ectopic_axon_AP_locs = axon_AP_locs_big{row,1};
    AP_data.evoked_axon_AP_locs{row,1} = evoked_axon_AP_locs(evoked_i);
    AP_data.ectopic_axon_AP_locs{row,1} = evoked_axon_AP_locs(ectopic_i);
    evoked_axon_AP_fidelity = axon_AP_fidelity_big{row,1};
    ectopic_axon_AP_fidelity = axon_AP_fidelity_big{row,1};
    AP_data.evoked_axon_AP_fidelity{row,1} = evoked_axon_AP_fidelity(evoked_i);
    AP_data.ectopic_axon_AP_fidelity{row,1} = evoked_axon_AP_fidelity(ectopic_i);
    evoked_axon_AP_latencies = axon_AP_latencies_big{row,1};
    ectopic_axon_AP_latencies = axon_AP_latencies_big{row,1};
    AP_data.evoked_axon_AP_latencies{row,1} = evoked_axon_AP_latencies(evoked_i);
    AP_data.ectopic_axon_AP_latencies{row,1} = evoked_axon_AP_latencies(ectopic_i);
    evoked_axon_AP_peaks = axon_AP_peaks_big{row,1};
    ectopic_axon_AP_peaks = axon_AP_peaks_big{row,1};
    AP_data.evoked_axon_AP_peaks{row,1} = evoked_axon_AP_peaks(evoked_i);
    AP_data.ectopic_axon_AP_peaks{row,1} = evoked_axon_AP_peaks(ectopic_i);

end

%% Find Indices of Axonal APs (Independently of Soma Peaks)

axon_AP_locs2 = []; % Create array to store indices of soma APs

for sweep = 1:size(data,3)

trace = data(:,axon_i_chan,sweep);

% Subtract resting membrane potential for peak finding
trace = -(trace - mean(trace(1:current_indices(1)-1)));

% Define region to calculate mean and SD of noise
noise_region = data(1:current_indices(1)-1,axon_i_chan,sweep);

%mean_noise = mean(noise_region); % get mean of noise
sd_noise = std(noise_region); % get SD of noise

% Set threhold for counting successful axon spike
axon_spike_threshold = 5*sd_noise;
% Set minimum time between peaks (i.e. minimum width of spike)
min_spike_width = 0.002; % 2 ms
% Set min on either side before the signal attains a higher value
axon_spike_prom = 10;

fs = 1e6 / si; % sampling rate
[amps,locs] = findpeaks(trace, fs, 'MinPeakHeight', axon_spike_threshold,...
    'MinPeakDistance', min_spike_width, 'MinPeakProminence', axon_spike_prom);
locs = int64(locs * fs);

% figure
% plot(trace)
% hold on
% %Plot triangles at the peak indices
% x = 1:height(trace);
% plot(x(locs), trace(locs), 'v', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

axon_AP_locs2 = [axon_AP_locs2; {locs}];

end

AP_data.axon_AP_locs2 = axon_AP_locs2; % Add to AP data table

current_injection_start = current_indices(1); % index where I injection start
current_injection_end = current_indices(1)+sweep_length*(1e6 / si); % index where I injection ends

% Sort APs into evoked and ectopic
for row = 1:height(AP_data)

    trace = AP_data.axon_AP_locs2{row,1};
    % get APs inside current injection
    evoked = trace(trace > current_injection_start & trace < current_injection_end);
    % get APs outside current injection
    EAPs = trace(trace < current_injection_start | trace > current_injection_end);
    
    AP_data.evoked_axon_AP_locs2{row,1} = evoked;
    AP_data.ectopic_axon_AP_locs2{row,1} = EAPs;

end

%% Count number of evoked and ectopic somatic and axonal APs per sweep

AP_data.num_soma_evoked_APs = cellfun(@numel, AP_data.evoked_soma_AP_locs);
AP_data.num_soma_EAPs = cellfun(@numel, AP_data.ectopic_soma_AP_locs);
AP_data.num_axon_evoked_APs = cellfun(@(x) sum(x == 1), AP_data.evoked_axon_AP_fidelity);
AP_data.num_axon_evoked_APs2 = cellfun(@numel, AP_data.evoked_axon_AP_locs2);
AP_data.num_axon_EAPs2 = cellfun(@numel, AP_data.ectopic_axon_AP_locs2);

