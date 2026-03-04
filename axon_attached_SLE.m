%% Axon attached signal analysis during bicuculline-induced seizure-like events
% Author: Sophie Liebergall
% Updated: 1/14/25
% Purpose: Extract signal from loose seal recordings of axonal APs during
% bicuculline-induced seizure-like events

%% Load in .abf file

%clear all

% Set path to .abf file
%pathFileName = "\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\data\axon_propagation\slice_electrophysiology\PV\bicuculline\for_analysis\243_2_SLE.abf";

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

trace = data(:,soma_v_chan,1); % Get soma voltage in trace

% Subtract resting membrane potential for peak finding
rest = mean(trace(1:1000));
trace = trace - rest;

% Display the electrophysiology trace
figure('WindowState', 'maximized'); % Create a figure maximized to full screen
plot(trace, 'b'); % Plot the trace in blue
xlabel('Sample Index');
ylabel('Amplitude');
title('Electrophysiology Trace');
hold on;

% Select a region to analyze
disp('Click and drag to select a region to analyze.');
h = drawrectangle; % Interactive rectangle selection
wait(h); % Wait for the user to complete the selection
region = round(h.Position); % Get the region coordinates as [x, y, width, height]
analyze_start = region(1); % Starting x-coordinate of the region
analyze_end = analyze_start + region(3); % Ending x-coordinate of the region
plot(analyze_start:analyze_end, trace(analyze_start:analyze_end), 'r'); % Highlight the selected region in red
disp(['Selected region: ', num2str(analyze_start), ' to ', num2str(analyze_end)]);

% Indicate the "SD_peak" point
% disp('Click to select the SD_peak point.');
% [x_SD_peak, ~] = ginput(1); % Select a single point
% SD_peak = round(x_SD_peak);
% plot(SD_peak, trace(SD_peak), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Mark the SD_peak point with a red circle
% disp(['SD_peak point: ', num2str(SD_peak)]);

% Hold off after plotting everything
hold off;

% Output the selected points and region
disp('Analysis Complete:');
disp(['Analyze Region: ', num2str(analyze_start), ' to ', num2str(analyze_end)]);
%disp(['SD_peak: ', num2str(SD_peak)]);

% Only look for spikes within set region
trace_crop = trace(analyze_start:analyze_end);

% Set threhold for counting successful soma spike
soma_spike_threshold = -20 - rest;
% Set minimum time between peaks (i.e. minimum width of spike)
min_spike_width = 0.002; % 2 ms
% Set minimum time between peaks (i.e. maximum width between spike)
max_spike_width = 0.010; % 2 ms
% Set min on either side before the signal attains a higher value
soma_spike_prom = 20;

fs = 1e6 / si; % sampling rate
[peaks,locs,widths,amps] = findpeaks(trace_crop, fs, 'MinPeakHeight', soma_spike_threshold,...
    'MinPeakDistance', min_spike_width, ...
    'WidthReference','halfprom','MaxPeakWidth',max_spike_width,...
    'MinPeakProminence', soma_spike_prom);
locs = int64(locs * fs);

% figure('WindowState', 'maximized'); % Create a figure maximized to full screen
% plot(trace_crop)
% hold on
% %Plot triangles at the peak indices
% x = (1:1:height(trace_crop))';
% plot(x(locs), trace_crop(locs), 'v', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
%xlim([2.39e7 2.4e7])
soma_AP_locs = locs + analyze_start;

AP_data.soma_AP_locs = soma_AP_locs; % Add to AP data table
AP_data.soma_AP_peaks = peaks;

close all;

%% Find Indices of Axonal APs
    
num_APs = size(soma_AP_locs,1); % number of APs in trace

% Define region to calculate mean and SD of noise
noise_region = data(1:1000,axon_i_chan,1);

%mean_noise = mean(noise_region); % get mean of noise
sd_noise = std(noise_region); % get SD of noise

% Set threhold for counting successful axon spike
axon_spike_threshold = 3*sd_noise;

% Create arrays to store data (each row is a sweep)
axon_AP_locs = [];
axon_AP_peaks = [];
axon_AP_latencies = [];
axon_AP_fidelity = [];

axon_AP_window = 0.75; % Set window (one sided) to detect axon apike (in ms)
axon_AP_window_i = axon_AP_window/(1e-3*si); % convert to indices

for spike = 1:size(soma_AP_locs,1) % loop over each spike in trace

    soma_AP_peak_i = soma_AP_locs(spike); % get index of soma AP peak
    
    % Select region to search for axon spike
    search_region = data(soma_AP_peak_i-axon_AP_window_i:soma_AP_peak_i+axon_AP_window_i,...
        axon_i_chan,1);

    % Get local minimum of search region
    [axon_AP_inward, axon_AP_peak_i] = min(search_region);
    axon_AP_peak_i = axon_AP_peak_i + (soma_AP_peak_i-axon_AP_window_i);
    % Get mean value 1 to 0.5 ms before AP peak
    baseline = mean(data(axon_AP_peak_i-(axon_AP_window_i*2):axon_AP_peak_i-(axon_AP_window_i),...
        axon_i_chan,1));
    axon_AP_peak = abs(axon_AP_inward - baseline);

    % Get first local max of outward current after axon spike peak
    width_window = 0.8; % Set width of window to look for max outward current after peak (in ms)
    width_window = width_window/(1e-3*si); % convert to indices
    width_search_region = data(axon_AP_peak_i:axon_AP_peak_i+width_window,axon_i_chan,1);
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
 
AP_data.axon_AP_locs = axon_AP_locs';
AP_data.axon_AP_fidelity = axon_AP_fidelity';
AP_data.axon_AP_latencies = axon_AP_latencies';
AP_data.axon_AP_peaks = axon_AP_peaks';

%%

% axon_trace = data(:,axon_i_chan,1);
% figure('WindowState', 'maximized'); % Create a figure maximized to full screen
% plot(axon_trace)
% hold on
% %Plot triangles at the peak indices
% x = (1:1:height(axon_trace))';
% plot(x(axon_AP_locs), axon_trace(axon_AP_locs), 'v', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% plot(trace)
% plot(x(soma_AP_locs), trace(soma_AP_locs), 'v', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
% xlim([2.4e7 2.405e7])

%% Divide the spikes into individual bursts

% Parameters
ISI_threshold_ms = 500; % Maximum inter-spike interval (in milliseconds) to consider spikes as part of the same burst

% Convert ISI threshold from ms to indices
ISI_threshold = ISI_threshold_ms / 1000 * fs;

% Find inter-spike intervals
spike_times = soma_AP_locs(:); % Ensure soma_AP_locs is a column vector
ISIs = diff(spike_times); % Calculate inter-spike intervals (in indices)

% Identify bursts
burst_indices = [1; find(ISIs > ISI_threshold) + 1; length(spike_times) + 1];
num_bursts = length(burst_indices) - 1;

% Assign bursts to spikes
burst_assignments = zeros(length(spike_times), 1); % Initialize with zeros
for i = 1:num_bursts
    burst_assignments(burst_indices(i):(burst_indices(i+1)-1)) = i; % Assign burst number
end

% % Optional: Plot bursts on the trace
% figure('WindowState', 'maximized');
% plot(trace, 'b'); % Plot the trace in blue
% hold on;
% for i = 1:num_bursts
%     burst_spikes = bursts{i};
%     plot(burst_spikes, trace(burst_spikes), 'o', 'MarkerSize', 8, 'DisplayName', ['Burst ' num2str(i)]);
% end
% xlabel('Sample Index');
% ylabel('Amplitude');
% title('Bursts on Electrophysiology Trace');
% legend('show');
% hold off;

AP_data.bursts = burst_assignments;
