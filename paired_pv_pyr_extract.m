%%% Extract Signal from PV-Pyr Paired Recording

% Author: Sophie Liebergall
% Updated: 9/10/25
% Purpose: Extract signal from paired recording of a PV-IN and Pyr cell

%% Load in .abf file

%clear all

% Set path to .abf file
%pathFileName = "\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\data\axon_propagation\slice_electrophysiology\PV_pyr_pairs\380_2_100Hz.abf";

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

pre_v_chan = 1;
post_i_chan = 2;
pre_i_chan = 3;
    
% Create vector with all indices in file
indices = 1:i;

% Create a vector with all indices converted to time points
time(indices,1) = (indices./(1e6./si));

% Create table to store somatic and axon AP data
AP_data = table();

%% Find Indices of Somatic APs

soma_AP_locs = []; % Create array to store indices of soma APs
soma_AP_peaks = []; % Create array to store soma AP peaks

for sweep = 1:size(data,3)

trace = data(:,pre_v_chan,sweep);

% Calculate resting membrane potential
rest = mean(trace(1:5000));

% Subtract resting membrane potential for peak finding
trace = trace - rest;

% Set threhold for counting successful soma spike
soma_spike_threshold = -10 - rest;
% Set minimum time between peaks (i.e. minimum width of spike)
min_spike_width = 0.002; % 2 ms
% Set maximum time between peaks (i.e. maximum width between spike)
max_spike_width = 0.010; % 2 ms
% Set min on either side before the signal attains a higher value
soma_spike_prom = 20;

fs = 1e6 / si; % sampling rate
[peaks,locs,widths,amps] = findpeaks(trace, fs, 'MinPeakHeight', soma_spike_threshold,...
    'MinPeakDistance', min_spike_width, ...
    'WidthReference','halfprom','MaxPeakWidth',max_spike_width,...
    'MinPeakProminence', soma_spike_prom);
locs = int64(locs * fs);

soma_AP_locs = [soma_AP_locs; {locs}];
soma_AP_peaks = [soma_AP_peaks; {peaks}];

end

AP_data.soma_AP_locs = soma_AP_locs; % Add to AP data table
AP_data.soma_AP_peaks = soma_AP_peaks;

%% Find Indices of Inward IPSCs

% Create arrays to store data across sweeps
 IPSC_locs_big = {};
 IPSC_fidelity_big = {};
 IPSC_latencies_big = {};
 IPSC_peaks_big = {};

for sweep = 1:size(data,3) % loop over each sweep in recording

soma_AP_locs = AP_data.soma_AP_locs{sweep,1};
    
num_APs = size(soma_AP_locs,1); % number of APs in trace

% Define region to calculate mean and SD of noise
noise_region = data(1:1000,post_i_chan,sweep);

%mean_noise = mean(noise_region); % get mean of noise
sd_noise = std(noise_region); % get SD of noise

% Set threhold for counting successful axon spike
IPSC_threshold = 5*sd_noise;

% Create arrays to store data (each row is a sweep)
IPSC_locs = [];
IPSC_peaks = [];
IPSC_latencies = [];
IPSC_fidelity = [];

IPSC_window = 2; % Set window (one sided) to detect IPSCs (in ms)
IPSC_window_i = IPSC_window/(1e-3*si); % convert to indices

% Get baseline for calculating IPSC peak amplitude
baseline = mean(noise_region);

for spike = 1:size(soma_AP_locs,1) % loop over each spike in trace

    soma_AP_peak_i = soma_AP_locs(spike); % get index of soma AP peak

    if soma_AP_peak_i < (height(data) - IPSC_window_i)
        IPSC_start_i = soma_AP_peak_i;
        IPSC_end_i = soma_AP_peak_i + IPSC_window_i;
    else
    % Select region to search for axon spike
        IPSC_start_i = soma_AP_peak_i;
        IPSC_end_i = height(data);
    end

    % Select region to search for axon spike
    search_region = data(IPSC_start_i:IPSC_end_i,...
        post_i_chan,sweep);
    % Get local minimum of search region
    [IPSC_inward, IPSC_peak_i] = min(search_region);
    

    if (soma_AP_peak_i > IPSC_window_i*2) && (soma_AP_peak_i < (height(data) - IPSC_window_i))
        IPSC_peak_i = IPSC_peak_i + (soma_AP_peak_i-IPSC_window_i);
        axon_mean_start_i = IPSC_peak_i-(IPSC_window_i*2);
        axon_mean_end_i = IPSC_peak_i-(IPSC_window_i);
    elseif (soma_AP_peak_i < IPSC_window_i)
        axon_mean_start_i = 1;
        axon_mean_end_i = 1+IPSC_window_i;
    else
        IPSC_peak_i = IPSC_peak_i + (soma_AP_peak_i-IPSC_window_i);
        axon_mean_start_i = IPSC_peak_i-(IPSC_window_i*2);
        axon_mean_end_i = IPSC_peak_i-(IPSC_window_i);
    end

    % Get mean value 1 to 0.5 ms before AP peak
    IPSC_peak = abs(IPSC_inward - baseline);
    
    IPSC_locs = [IPSC_locs IPSC_peak_i];
    IPSC_peaks = [IPSC_peaks IPSC_peak];

    if IPSC_peak > IPSC_threshold
        IPSC_fidelity = [IPSC_fidelity 1];
    else
        IPSC_fidelity = [IPSC_fidelity 0];
    end

    % Calculate latency from somatic spike peak to axon spike peak
    latency = IPSC_peak_i - soma_AP_peak_i;
    latency = double(latency) * (1e-3*si); % convert to ms
    IPSC_latencies = [IPSC_latencies latency];

end
 
 IPSC_locs_big{sweep, 1} = IPSC_locs;
 IPSC_fidelity_big{sweep, 1} = IPSC_fidelity;
 IPSC_latencies_big{sweep, 1} = IPSC_latencies;
 IPSC_peaks_big{sweep, 1} = IPSC_peaks;
 

end

AP_data.IPSC_locs = IPSC_locs_big;
AP_data.IPSC_fidelity = IPSC_fidelity_big;
AP_data.IPSC_latencies = IPSC_latencies_big;
AP_data.IPSC_peaks = IPSC_peaks_big;

%% Count number of evoked and ectopic somatic and axonal APs per sweep

AP_data.num_soma_APs = cellfun(@numel, AP_data.soma_AP_locs);
AP_data.num_IPSCs = cellfun(@(x) sum(x == 1), AP_data.IPSC_fidelity);


