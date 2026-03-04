%%% Model Rundown in Axonal Action Potential Fidelity

% Author: Sophie Liebergall
% Updated: 12/5/24
% Purpose: Analyze "rundown" of conduction during normal IC steps protocols

% Inputs: AP_data table (output from axon_attached_IC_steps_EAPs)

clearvars %-except AP_data

% Load in axon AP data table for PV-INs
% IC steps data
ICsteps_data = load("\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\analysis\axon_propagation\axonal_fidelity\PV_all_axon_AP_data.mat");
ICsteps_data = ICsteps_data.all_AP_data;

% IC pulses data
pulses_data = load("\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\analysis\axon_propagation\axonal_fidelity\PV_all_axon_pulse_AP_data.mat");
pulses_data = pulses_data.all_AP_data;

% IC pulses data
var_pulses_data = load("\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\analysis\axon_propagation\axonal_fidelity\PV_all_axon_var_pulse_AP_data.mat");
var_pulses_data = var_pulses_data.all_AP_data;

% Read in intrinsic properties and morphology data
morph_data_path = "\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\analysis\axon_propagation\PV_intrinsic_props_morphology_V1.xlsx";
morph_data = readtable(morph_data_path);

% si = sampling interval in microseconds (us)
%       Note: sampling rate (Hz) = 1e6 / si (us)
si = 10;
fs = 1e6 ./ si; % sampling rate (Hz)
sweep_length = 3; % sweep length (sec)

% Set default text properties to use Arial font
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultUicontrolFontName', 'Arial');

% Set colors
soma_evoked_color = [128/255 130/255 133/255];
axon_evoked_color = '#bcbec0';
soma_ectopic_color = [185/255 37/255 121/255];
axon_ectopic_color = [230/255 108/255 177/255];
ci_color = [252/255 222/255 156/225];
pv_color = [68/255 119/255 170/225];
ndnf_color = [34/255 136/255 51/225];

% Set directoy to save figures
saveDir = "\\ressmb01.research.chop.edu\Goldberg_Lab\Goldberg Lab - Personal Folders\Sophie Liebergall\analysis\axon_propagation\axonal_fidelity\PV_ICsteps_fidelity_figs\decay_model";
%fileName = "7392_1";

%% Create table to store data with all sweeps combined

% Table to store data for all cells for IC steps protocol
ICsteps_T = table();

for i = 1:size(ICsteps_data,2)
        
% Get AP data for individual cell
AP_data = ICsteps_data(i).table;

% Create table to store data for cell
T = table();

all_soma_AP_locs = []; % array to store all soma AP locs
all_axon_AP_peaks = []; % array to store all axon AP peaks
all_latencies = []; % array to store soma-axon AP latencies

% Sweep length in locs
sweep_length_i = sweep_length * fs;

for row = 1:height(AP_data) % loop over each sweep

    sweep_locs = cell2mat(AP_data.soma_AP_locs(row));
    sweep_locs = sweep_locs + (sweep_length_i*(row-1)); % add indices from previous sweeps

    all_soma_AP_locs = [all_soma_AP_locs; sweep_locs]; % accumulate in a table
    all_axon_AP_peaks = [all_axon_AP_peaks; cell2mat(AP_data.axon_AP_peaks(row))'];
    all_latencies = [all_latencies; cell2mat(AP_data.axon_AP_latencies(row))'];

end

% Convert indices to time (seconds)
all_soma_AP_locs = double(all_soma_AP_locs) ./ fs;

T.soma_AP_locs = all_soma_AP_locs;
T.axon_AP_peaks = all_axon_AP_peaks;
T.axon_AP_peaks_norm = all_axon_AP_peaks./max(all_axon_AP_peaks);
T.latencies = all_latencies;
T.mouse = repmat(double(ICsteps_data(i).mouse), height(all_soma_AP_locs), 1);
T.cell = repmat(double(ICsteps_data(i).cell), height(all_soma_AP_locs), 1);
T.ID = repmat(string(ICsteps_data(i).mouse)+"_"+string(ICsteps_data(i).cell), height(all_soma_AP_locs), 1);

ICsteps_T = [ICsteps_T; T];

end

% Table to store data for all cells for IC pulses protocol
pulses_T = table();

for i = 1:size(pulses_data,2)
        
% Get AP data for individual cell
AP_data = pulses_data(i).table;

% Create table to store data for cell
T = table();

all_soma_AP_locs = []; % array to store all soma AP locs
all_axon_AP_peaks = []; % array to store all axon AP peaks
all_latencies = []; % array to store soma-axon AP latencies

% Sweep length in locs
sweep_length_i = sweep_length * fs;

for row = 1:height(AP_data) % loop over each sweep

    sweep_locs = cell2mat(AP_data.soma_AP_locs(row));
    sweep_locs = sweep_locs + (sweep_length_i*(row-1)); % add indices from previous sweeps

    all_soma_AP_locs = [all_soma_AP_locs; sweep_locs]; % accumulate in a table
    all_axon_AP_peaks = [all_axon_AP_peaks; cell2mat(AP_data.axon_AP_peaks(row))'];
    all_latencies = [all_latencies; cell2mat(AP_data.axon_AP_latencies(row))'];

end

% Convert indices to time (seconds)
all_soma_AP_locs = double(all_soma_AP_locs) ./ fs;

T.soma_AP_locs = all_soma_AP_locs;
T.axon_AP_peaks = all_axon_AP_peaks;
T.axon_AP_peaks_norm = all_axon_AP_peaks./max(all_axon_AP_peaks);
T.latencies = all_latencies;
T.mouse = repmat(double(pulses_data(i).mouse), height(all_soma_AP_locs), 1);
T.cell = repmat(double(pulses_data(i).cell), height(all_soma_AP_locs), 1);
T.frequency = repmat(double(pulses_data(i).frequency), height(all_soma_AP_locs), 1);
T.ID = repmat(string(pulses_data(i).mouse)+"_"+string(pulses_data(i).cell)+"_"+string(pulses_data(i).frequency), height(all_soma_AP_locs), 1);

pulses_T = [pulses_T; T];

end

% Table to store data for all cells for IC pulses protocol
var_pulses_T = table();

for i = 1:size(var_pulses_data,2)
        
% Get AP data for individual cell
AP_data = var_pulses_data(i).table;

% Create table to store data for cell
T = table();

all_soma_AP_locs = []; % array to store all soma AP locs
all_axon_AP_peaks = []; % array to store all axon AP peaks
all_latencies = []; % array to store soma-axon AP latencies

% Sweep length in locs
sweep_length_i = sweep_length * fs;

for row = 1:height(AP_data) % loop over each sweep

    sweep_locs = cell2mat(AP_data.soma_AP_locs(row));
    sweep_locs = sweep_locs + (sweep_length_i*(row-1)); % add indices from previous sweeps

    all_soma_AP_locs = [all_soma_AP_locs; sweep_locs]; % accumulate in a table
    all_axon_AP_peaks = [all_axon_AP_peaks; cell2mat(AP_data.axon_AP_peaks(row))'];
    all_latencies = [all_latencies; cell2mat(AP_data.axon_AP_latencies(row))'];

end

% Convert indices to time (seconds)
all_soma_AP_locs = double(all_soma_AP_locs) ./ fs;

T.soma_AP_locs = all_soma_AP_locs;
T.axon_AP_peaks = all_axon_AP_peaks;
T.latencies = all_latencies;
T.axon_AP_peaks_norm = all_axon_AP_peaks./max(all_axon_AP_peaks);
T.mouse = repmat(double(var_pulses_data(i).mouse), height(all_soma_AP_locs), 1);
T.cell = repmat(double(var_pulses_data(i).cell), height(all_soma_AP_locs), 1);
T.ID = repmat(string(var_pulses_data(i).mouse)+"_"+string(var_pulses_data(i).cell), height(all_soma_AP_locs), 1);

var_pulses_T = [var_pulses_T; T];

end

% %% Model for Activity-dependent Decreases in Axonal Fidelity (IC steps)
% 
% % Load your data
% % Assume you have:
% % time_AP: Nx1 vector of AP times (in ms)
% % amp_axon: Nx1 vector of recorded axonal AP amplitudes
% % amp_soma: Nx1 vector of corresponding soma AP amplitudes
% % N_APs: Nx1 vector of cumulative AP count at each time point
% 
% ID_i = "8305_1"; % Select recording for modeling data
% 
% time_AP = ICsteps_T.soma_AP_locs(ICsteps_T.ID == ID_i); % time of APs in seconds
% amp_axon = ICsteps_T.axon_AP_peaks_norm(ICsteps_T.ID == ID_i); % vector of normalized axon AP peaks
% 
% % Define the model function
% decay_model = @(params, t, t_APs) params(1) * sum(exp(-(t - t_APs) / params(2)));
% 
% % Define the objective function for fitting
% objective_fun = @(params) amp_axon - arrayfun(@(t) ...
%     sum(decay_model(params, t, time_AP(time_AP < t))), time_AP, 'UniformOutput', true);
% 
% params0 = [max(amp_axon), 10]; % Initial guesses: A0 = max amplitude, tau = 10ms
% lb = [0, 0]; % Lower bounds (A0 and tau must be positive)
% ub = [Inf, Inf]; % Upper bounds
% 
% params_fit = lsqnonlin(objective_fun, params0, lb, ub);
% 
% % Extract fitted parameters
% A0_fit = params_fit(1);
% tau_fit = params_fit(2);
% 
% amp_axon_pred = arrayfun(@(t) ...
%     sum(decay_model(params_fit, t, time_AP(time_AP < t))), time_AP, 'UniformOutput', true);
% 
% figure;
% scatter(time_AP, amp_axon, 'bo'); hold on; % Original data
% plot(time_AP, amp_axon_pred, 'r-', 'LineWidth', 2); % Fitted model
% xlabel('Time (ms)');
% ylabel('AP Amplitude');
% legend('Data', 'Fitted Model');
% title('Axonal AP Amplitude Decay Model');
% hold off;

%% Decay Model for Activity-dependent Decreases in Axonal Fidelity (IC Steps)

num_cells = height(unique(ICsteps_T.ID));

params_fit_all = zeros(num_cells, 3); % Store fitted parameters for each cell
R2_all = zeros(num_cells, 1); % Store R² values

IDs = unique(ICsteps_T.ID); % array of file IDs

for i = 1:num_cells

ID_i = IDs(i); % Get recording ID of interest

time_AP = ICsteps_T.soma_AP_locs(ICsteps_T.ID == ID_i); % time of APs in seconds
amp_axon = ICsteps_T.axon_AP_peaks_norm(ICsteps_T.ID == ID_i); % vector of normalized axon AP peaks

% Define the corrected decay model
decay_model = @(params, t, t_APs) params(1) - params(2) * sum(exp(-(t - t_APs) / params(3)));

% Initial guesses for parameters [A0, Delta_A, tau]
params0 = [max(amp_axon), max(amp_axon) / 5, 10]; % Adjusted initial guesses

% Define the objective function for fitting
objective_fun = @(params) amp_axon - arrayfun(@(t) ...
    decay_model(params, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);

% Fit using nonlinear least squares
params_fit = lsqnonlin(objective_fun, params0, [0, 0, 0], [Inf, 100, 200]);
params_fit_all(i, :) = params_fit;

% Extract fitted parameters
A0_fit = params_fit(1);
Delta_A_fit = params_fit(2);
tau_fit = params_fit(3);

% Predict axonal amplitude using the fitted model
amp_axon_pred = arrayfun(@(t) ...
    decay_model(params_fit, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);

% Plot results
figure;
scatter(time_AP, amp_axon, 50, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color,'LineWidth',0.25,...
    'MarkerFaceAlpha',0.25); hold on; % Original data
plot(time_AP, amp_axon_pred, 'k-', 'LineWidth', 2); % Fitted model
xlabel('Time (ms)');
ylabel('Normalized axon AP amplitude');
legend('Data', 'Fitted Model');
%title('Axonal AP Amplitude Decay Model');
hold off;

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "decay_model_fit_" + ID_i;

%saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
%saveas(gcf, fullfile(saveDir,(tit)), 'png')

close all;

% Compute R²
    SS_res = sum((amp_axon - amp_axon_pred).^2);
    SS_tot = sum((amp_axon - mean(amp_axon)).^2);
    R2_all(i) = 1 - (SS_res / SS_tot);

end

% %% Decay Model 2 for Activity-dependent Decreases in Axonal Fidelity (IC Steps)
% 
% num_cells = height(unique(ICsteps_T.ID));
% 
% params_fit_all2 = zeros(num_cells, 2); % Store fitted parameters for each cell
% R2_all2 = zeros(num_cells, 1); % Store R² values
% 
% IDs = unique(ICsteps_T.ID); % array of file IDs
% 
% for i = 1:num_cells
% 
% ID_i = IDs(i); % Get recording ID of interest
% 
% time_AP = ICsteps_T.soma_AP_locs(ICsteps_T.ID == ID_i); % time of APs in seconds
% amp_axon = ICsteps_T.axon_AP_peaks_norm(ICsteps_T.ID == ID_i); % vector of normalized axon AP peaks
% amp0 = amp_axon(1); % first AP amplitude
% 
% % Define the corrected multiplicative decay model
% decay_model = @(params, t, t_APs, amp0) ...
%     amp0 * (isempty(t_APs) + ~isempty(t_APs) * prod(exp(-(t - t_APs) / params(2))));
% 
% % Initial guess for parameters [A0, tau]
% params0 = [max(amp_axon), 10];  % You may adjust based on your data
% 
% % Define the objective function for fitting
% objective_fun = @(params) amp_axon - arrayfun(@(t) ...
%     decay_model(params, t, time_AP(time_AP < t), amp0), time_AP, 'UniformOutput', true);
% 
% % Fit using nonlinear least squares
% params_fit = lsqnonlin(objective_fun, params0, [0, 0], [1, Inf]);
% params_fit_all2(i, :) = params_fit;
% 
% % Extract the fitted parameters
% A0_fit = params_fit(1);
% tau_fit = params_fit(2);
% 
% % Predict the amplitudes using the fitted model
% amp_axon_pred = arrayfun(@(t) decay_model(params_fit, t, time_AP(time_AP < t), amp0), time_AP);
% 
% % Plot results
% figure;
% scatter(time_AP, amp_axon, 50, 'MarkerFaceColor', pv_color,...
%     'MarkerEdgeColor', pv_color,'LineWidth',0.25,...
%     'MarkerFaceAlpha',0.25); hold on; % Original data
% plot(time_AP, amp_axon_pred, 'k-', 'LineWidth', 2); % Fitted model
% xlabel('Time (ms)');
% ylabel('Normalized axon AP amplitude');
% legend('Data', 'Fitted Model');
% %title('Axonal AP Amplitude Decay Model');
% hold off;
% 
% % Set the figure size
% set(gcf, 'Position', [100, 100, 800, 600]);
% 
% tit = "decay_model2_fit_" + ID_i;
% 
% saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
% saveas(gcf, fullfile(saveDir,(tit)), 'png')
% 
% close all;
% 
% % Compute R²
%     SS_res = sum((amp_axon - amp_axon_pred).^2);
%     SS_tot = sum((amp_axon - mean(amp_axon)).^2);
%     R2_all(i) = 1 - (SS_res / SS_tot);
% 
% end

%% Get general fitting parameters for IC steps and plot

med_params = median(params_fit_all, 1); % Median of all parameters

R2_all_med = zeros(num_cells, 1); % Store R² values

IDs = unique(ICsteps_T.ID); % array of file IDs

for i = 1:num_cells

ID_i = IDs(i); % Get recording ID of interest

time_AP = ICsteps_T.soma_AP_locs(ICsteps_T.ID == ID_i); % time of APs in seconds
amp_axon = ICsteps_T.axon_AP_peaks_norm(ICsteps_T.ID == ID_i); % vector of normalized axon AP peaks

% Define the corrected decay model
decay_model = @(params, t, t_APs) params(1) - params(2) * sum(exp(-(t - t_APs) / params(3)));

% Extract fitted parameters
A0_fit = params_fit(1);
Delta_A_fit = params_fit(2);
tau_fit = params_fit(3);

% Predict axonal amplitude using the fitted model
amp_axon_pred = arrayfun(@(t) ...
    decay_model(med_params, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);

% Plot results
figure;
scatter(time_AP, amp_axon, 50, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color,'LineWidth',0.25,...
    'MarkerFaceAlpha',0.25); hold on; % Original data
plot(time_AP, amp_axon_pred, 'k-', 'LineWidth', 2); % Fitted model
xlabel('Time (ms)');
ylabel('Normalized axon AP amplitude');
legend('Data', 'Fitted Model');
%title('Axonal AP Amplitude Decay Model');
hold off;

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "decay_model_fit_" + ID_i + "_med_params";

%saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
%saveas(gcf, fullfile(saveDir,(tit)), 'png')

close all;

% Compute R²
    SS_res = sum((amp_axon - amp_axon_pred).^2);
    SS_tot = sum((amp_axon - mean(amp_axon)).^2);
    R2_all_med(i) = 1 - (SS_res / SS_tot);

end

%% Decay Model for Activity-dependent Decreases in Axonal Fidelity (Pulses)

num_cells = height(unique(pulses_T.ID));

pulses_params_fit_all = zeros(num_cells, 3); % Store fitted parameters for each cell
pulses_R2_all = zeros(num_cells, 1); % Store R² values

IDs = unique(pulses_T.ID); % array of file IDs

for i = 1:num_cells

ID_i = IDs(i); % Get recording ID of interest

time_AP = pulses_T.soma_AP_locs(pulses_T.ID == ID_i); % time of APs in seconds
amp_axon = pulses_T.axon_AP_peaks_norm(pulses_T.ID == ID_i); % vector of normalized axon AP peaks

% Define the corrected decay model
decay_model = @(params, t, t_APs) params(1) - params(2) * sum(exp(-(t - t_APs) / params(3)));

% Initial guesses for parameters [A0, Delta_A, tau]
params0 = [max(amp_axon), max(amp_axon) / 5, 10]; % Adjusted initial guesses

% Define the objective function for fitting
objective_fun = @(params) amp_axon - arrayfun(@(t) ...
    decay_model(params, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);

% Fit using nonlinear least squares
params_fit = lsqnonlin(objective_fun, params0, [0, 0, 0], [Inf, 100, 200]);
pulses_params_fit_all(i, :) = params_fit;

% Extract fitted parameters
A0_fit = params_fit(1);
Delta_A_fit = params_fit(2);
tau_fit = params_fit(3);

% Predict axonal amplitude using the fitted model
amp_axon_pred = arrayfun(@(t) ...
    decay_model(params_fit, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);

% Plot results
figure;
scatter(time_AP, amp_axon, 50, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color,'LineWidth',0.25,...
    'MarkerFaceAlpha',0.25); hold on; % Original data
plot(time_AP, amp_axon_pred, 'k-', 'LineWidth', 2); % Fitted model
xlabel('Time (ms)');
ylabel('Normalized axon AP amplitude');
legend('Data', 'Fitted Model');
%title('Axonal AP Amplitude Decay Model');
hold off;

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "decay_model_fit_" + ID_i + "_pulses";

%saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
%saveas(gcf, fullfile(saveDir,(tit)), 'png')

close all;

% Compute R²
    SS_res = sum((amp_axon - amp_axon_pred).^2);
    SS_tot = sum((amp_axon - mean(amp_axon)).^2);
    pulses_R2_all(i) = 1 - (SS_res / SS_tot);

end

%% Decay Model for Activity-dependent Decreases in Axonal Fidelity (Variable Pulses)

num_cells = height(unique(var_pulses_T.ID));

var_pulses_params_fit_all = zeros(num_cells, 3); % Store fitted parameters for each cell
var_pulses_R2_all = zeros(num_cells, 1); % Store R² values

IDs = unique(var_pulses_T.ID); % array of file IDs

for i = 1:num_cells

ID_i = IDs(i); % Get recording ID of interest

time_AP = var_pulses_T.soma_AP_locs(var_pulses_T.ID == ID_i); % time of APs in seconds
amp_axon = var_pulses_T.axon_AP_peaks_norm(var_pulses_T.ID == ID_i); % vector of normalized axon AP peaks

% Define the corrected decay model
decay_model = @(params, t, t_APs) params(1) - params(2) * sum(exp(-(t - t_APs) / params(3)));

% Initial guesses for parameters [A0, Delta_A, tau]
params0 = [max(amp_axon), max(amp_axon) / 5, 10]; % Adjusted initial guesses

% Define the objective function for fitting
objective_fun = @(params) amp_axon - arrayfun(@(t) ...
    decay_model(params, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);

% Fit using nonlinear least squares
params_fit = lsqnonlin(objective_fun, params0, [0, 0, 0], [Inf, 100, 50]);
var_pulses_params_fit_all(i, :) = params_fit;

% Extract fitted parameters
A0_fit = params_fit(1);
Delta_A_fit = params_fit(2);
tau_fit = params_fit(3);

% Predict axonal amplitude using the fitted model
amp_axon_pred = arrayfun(@(t) ...
    decay_model(params_fit, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);

% Plot results
figure;
scatter(time_AP, amp_axon, 50, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color,'LineWidth',0.25,...
    'MarkerFaceAlpha',0.25); hold on; % Original data
plot(time_AP, amp_axon_pred, 'k-', 'LineWidth', 2); % Fitted model
xlabel('Time (ms)');
ylabel('Normalized axon AP amplitude');
legend('Data', 'Fitted Model');
%title('Axonal AP Amplitude Decay Model');
hold off;

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "decay_model_fit_" + ID_i + "_var_pulses";

%saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
%saveas(gcf, fullfile(saveDir,(tit)), 'png')

close all;

% Compute R²
    SS_res = sum((amp_axon - amp_axon_pred).^2);
    SS_tot = sum((amp_axon - mean(amp_axon)).^2);
    var_pulses_R2_all(i) = 1 - (SS_res / SS_tot);

end

% %% Model for Activity-dependent Decreases in Axonal Fidelity (Raw Signal)
% 
% num_cells = height(unique(ICsteps_T.ID));
% 
% params_fit_all2 = zeros(num_cells, 3); % Store fitted parameters for each cell
% R2_all2 = zeros(num_cells, 1); % Store R² values
% 
% IDs = unique(ICsteps_T.ID); % array of file IDs
% 
% for i = 1:num_cells
% 
% ID_i = IDs(i); % Get recording ID of interest
% 
% time_AP = ICsteps_T.soma_AP_locs(ICsteps_T.ID == ID_i); % time of APs in seconds
% amp_axon = ICsteps_T.axon_AP_peaks(ICsteps_T.ID == ID_i); % vector of normalized axon AP peaks
% 
% % Define the corrected decay model
% decay_model = @(params, t, t_APs) params(1) - params(2) * sum(exp(-(t - t_APs) / params(3)));
% 
% % Initial guesses for parameters [A0, Delta_A, tau]
% params0 = [max(amp_axon), max(amp_axon) / 5, 10]; % Adjusted initial guesses
% 
% % Define the objective function for fitting
% objective_fun = @(params) amp_axon - arrayfun(@(t) ...
%     decay_model(params, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);
% 
% % Fit using nonlinear least squares
% params_fit = lsqnonlin(objective_fun, params0, [0, 0, 0], [Inf, Inf, Inf]);
% params_fit_all(i, :) = params_fit;
% 
% % Extract fitted parameters
% A0_fit = params_fit(1);
% Delta_A_fit = params_fit(2);
% tau_fit = params_fit(3);
% 
% % Predict axonal amplitude using the fitted model
% amp_axon_pred = arrayfun(@(t) ...
%     decay_model(params_fit, t, time_AP(time_AP < t)), time_AP, 'UniformOutput', true);
% 
% % Plot results
% figure;
% scatter(time_AP, amp_axon, 50, 'MarkerFaceColor', pv_color,...
%     'MarkerEdgeColor', pv_color,'LineWidth',0.25,...
%     'MarkerFaceAlpha',0.25); hold on; % Original data
% plot(time_AP, amp_axon_pred, 'k-', 'LineWidth', 2); % Fitted model
% xlabel('Time (ms)');
% ylabel('Normalized axon AP amplitude');
% legend('Data', 'Fitted Model');
% %title('Axonal AP Amplitude Decay Model');
% hold off;
% 
% % Set the figure size
% set(gcf, 'Position', [100, 100, 800, 600]);
% 
% tit = "decay_model_fit_raw_" + ID_i;
% 
% saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
% saveas(gcf, fullfile(saveDir,(tit)), 'png')
% 
% close all;
% 
% % Compute R²
%     SS_res = sum((amp_axon - amp_axon_pred).^2);
%     SS_tot = sum((amp_axon - mean(amp_axon)).^2);
%     R2_all2(i) = 1 - (SS_res / SS_tot);
% 
% end

%% Correlate fit parameters with cell properties

%%% Create table with params from fitting (IC steps and pulse protocols)

% Create table of IC steps data for merging
ICsteps_comb_T = table();
ICsteps_comb_T.ID = unique(ICsteps_T.ID);
ICsteps_comb_T.param1 = params_fit_all(:,1);
ICsteps_comb_T.param2 = params_fit_all(:,2);
ICsteps_comb_T.param3 = params_fit_all(:,3);
ICsteps_comb_T.R2 = R2_all;
ICsteps_comb_T.protocol = repmat("ICsteps", height(ICsteps_comb_T), 1);
max_axon_AP = varfun(@max, ICsteps_T, 'GroupingVariables', 'ID', 'InputVariables', 'axon_AP_peaks');
ICsteps_comb_T.max_axon_AP = max_axon_AP.max_axon_AP_peaks;
latencies = varfun(@mean, ICsteps_T, 'GroupingVariables', 'ID', 'InputVariables', 'latencies');
ICsteps_comb_T.AP_latency = latencies.mean_latencies;

% Create table of pulses data for merging
pulses_comb_T = table();
pulses_comb_T.ID_freq = unique(pulses_T.ID);
pulses_comb_T.param1 = pulses_params_fit_all(:,1);
pulses_comb_T.param2 = pulses_params_fit_all(:,2);
pulses_comb_T.param3 = pulses_params_fit_all(:,3);
pulses_comb_T.R2 = pulses_R2_all;
pulses_comb_T.ID = split(pulses_comb_T.ID_freq, '_'); % Split at '&'
pulses_comb_T.ID = pulses_comb_T.ID(:,1) + "_" + pulses_comb_T.ID(:,2); % Keep only the first part
pulses_comb_T.protocol = repmat("pulses", height(pulses_comb_T), 1);
max_axon_AP = varfun(@max, pulses_T, 'GroupingVariables', 'ID', 'InputVariables', 'axon_AP_peaks');
pulses_comb_T.max_axon_AP = max_axon_AP.max_axon_AP_peaks;
latencies = varfun(@mean, pulses_T, 'GroupingVariables', 'ID', 'InputVariables', 'latencies');
pulses_comb_T.AP_latency = latencies.mean_latencies;

% Create table of variable pulses data for merging
var_pulses_comb_T = table();
var_pulses_comb_T.ID_freq = unique(var_pulses_T.ID);
var_pulses_comb_T.param1 = var_pulses_params_fit_all(:,1);
var_pulses_comb_T.param2 = var_pulses_params_fit_all(:,2);
var_pulses_comb_T.param3 = var_pulses_params_fit_all(:,3);
var_pulses_comb_T.R2 = var_pulses_R2_all;
var_pulses_comb_T.ID = split(var_pulses_comb_T.ID_freq, '_'); % Split at '&'
var_pulses_comb_T.ID = var_pulses_comb_T.ID(:,1) + "_" + var_pulses_comb_T.ID(:,2); % Keep only the first part
var_pulses_comb_T.protocol = repmat("var_pulses", height(var_pulses_comb_T), 1);
max_axon_AP = varfun(@max, var_pulses_T, 'GroupingVariables', 'ID', 'InputVariables', 'axon_AP_peaks');
var_pulses_comb_T.max_axon_AP = max_axon_AP.max_axon_AP_peaks;
latencies = varfun(@mean, var_pulses_T, 'GroupingVariables', 'ID', 'InputVariables', 'latencies');
var_pulses_comb_T.AP_latency = latencies.mean_latencies;

% Create array with desired variables in corr table
selectedVars = {'ID', 'param1', 'param2','param3', 'R2', 'protocol', 'max_axon_AP', 'AP_latency'};

% Combine table into correlations table
comb_T = vertcat(ICsteps_comb_T(:, selectedVars),...
    pulses_comb_T(:, selectedVars), var_pulses_comb_T(:, selectedVars));

%%% Add ephys and morphology data for each cell
% Add ID column to morph_T
%morph_data.ID = string(morph_data.mouse) + "_" + string(morph_data.cell);

% Find matching rows in morph_data based on ID
%[found, idx] = ismember(comb_T.ID, morph_data.ID);

corr_T = comb_T;

% Loop through all variables in morph_data (except ID) and assign values
varNames = setdiff(morph_data.Properties.VariableNames, 'ID'); % Get variable names excluding ID

% for i = 1:length(varNames)
%     % Assign values only for matching rows, keeping the original order
%     corr_T.(varNames{i})(found) = morph_data.(varNames{i})(idx(found));
% end

%% Plot fitting parameters

% Example: Assume corr_T is your table
% Columns: ID, param1, param2, param3, R2, protocol

% Ensure 'protocol' is categorical (for coloring)
corr_T.protocol = categorical(corr_T.protocol);

% Get variable names of interest
vars = {'param1', 'param2', 'param3', 'R2'};
varNames = {'A0','deltaA','tau','R^2'};

% Create a figure with subplots for each parameter
% figure;
% tiledlayout(2,2, 'TileSpacing', 'compact');

% Define colors for each protocol
protocols = categories(corr_T.protocol);
numProtocols = numel(protocols);
colors = lines(numProtocols); % or use other colormaps, e.g., 'parula', 'turbo'

close all;

for i = 1:numel(vars)
    figure;
    hold on

    % Plot a global boxplot *at a separate x-position* (e.g., after protocols)
    x_box = 1;
    boxplot(corr_T{:, vars{i}}, 'Positions', x_box, 'Colors', 'k',...
        'Symbol', '','Widths', 0.4, 'MedianStyle', 'line');

    % Loop over protocols and plot each swarm
    for j = 1:numProtocols
        protocol_j = protocols{j};
        idx = corr_T.protocol == protocol_j;
        s = swarmchart(ones(sum(idx),1), corr_T{idx, vars{i}}, ...
            100, repmat(colors(j,:), sum(idx),1), 'filled', ...
            'MarkerFaceAlpha', 0.2,'MarkerEdgeColor','k',...
            'MarkerEdgeAlpha',0.3);
        s.XJitterWidth = 0.3;
    end
    
    ylabel(varNames{i});
    xticks([]); % No x-ticks for individual subplots
    title(vars{i});
    xlim([0.5 1.5])
    ylim([0 (1.25*max(corr_T{:, vars{i}}))])

    % Add a legend
    lgd = legend({'Steps', 'Regular pulses', 'Variable pulses'}, 'Location', 'bestoutside');
    title(lgd, 'Protocol');

    % Set the figure size
    set(gcf, 'Position', [100, 100, 400, 600]);

    tit = varNames{i} + "_values";

    saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
    saveas(gcf, fullfile(saveDir,(tit)), 'png')

end

%% Plot fitting parameters - tau outliers excluded

% Example: Assume corr_T is your table
% Columns: ID, param1, param2, param3, R2, protocol

corr_T.protocol = categorical(corr_T.protocol);

% Define colors for each protocol
protocols = categories(corr_T.protocol);
numProtocols = numel(protocols);
colors = lines(numProtocols); % or use other colormaps, e.g., 'parula', 'turbo'

corr_T_plot = corr_T(corr_T.param3 < 199,:);

close all;


    figure;
    hold on

    % Plot a global boxplot *at a separate x-position* (e.g., after protocols)
    x_box = 1;
    boxplot(corr_T_plot{:, 'param3'}, 'Positions', x_box, 'Colors', 'k',...
        'Symbol', '','Widths', 0.4, 'MedianStyle', 'line');

    % Loop over protocols and plot each swarm
    for j = 1:numProtocols
        protocol_j = protocols{j};
        idx = corr_T_plot.protocol == protocol_j;
        s = swarmchart(ones(sum(idx),1), corr_T_plot{idx, 'param3'}, ...
            100, repmat(colors(j,:), sum(idx),1), 'filled', ...
            'MarkerFaceAlpha', 0.2,'MarkerEdgeColor','k',...
            'MarkerEdgeAlpha',0.3);
        s.XJitterWidth = 0.3;
    end
    
    ylabel('tau');
    xticks([]); % No x-ticks for individual subplots
    title('tau');
    xlim([0.5 1.5])
    ylim([0 (1.25*max(corr_T_plot{:, 'param3'}))])

    % Add a legend
    lgd = legend({'Steps', 'Regular pulses', 'Variable pulses'}, 'Location', 'bestoutside');
    title(lgd, 'Protocol');

    % Set the figure size
    set(gcf, 'Position', [100, 100, 400, 600]);

    tit = "tau_boundsexcluded_values";

    saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
    saveas(gcf, fullfile(saveDir,(tit)), 'png')

    %% Plot fitting parameters - outliers excluded

% Example: Assume corr_T is your table
% Columns: ID, param1, param2, param3, R2, protocol

% Ensure 'protocol' is categorical (for coloring)
corr_T.protocol = categorical(corr_T.protocol);

% Get variable names of interest
vars = {'param1', 'param2', 'param3', 'R2'};
varNames = {'A0','deltaA','tau','R^2'};

% Create a figure with subplots for each parameter
% figure;
% tiledlayout(2,2, 'TileSpacing', 'compact');

% Define colors for each protocol
protocols = categories(corr_T.protocol);
numProtocols = numel(protocols);
colors = lines(numProtocols); % or use other colormaps, e.g., 'parula', 'turbo'

close all;

for i = 1:numel(vars)
    figure;
    hold on

    raw = corr_T{:, vars{i}};
    Q1 = quantile(raw, 0.25);
    Q3 = quantile(raw, 0.75);
    IQR = Q3 - Q1;

    plot_corr_T = corr_T(raw >= (Q1 - 1.5*IQR) & raw <= (Q3 + 1.5*IQR),:);

    % Plot a global boxplot *at a separate x-position* (e.g., after protocols)
    x_box = 1;
    boxplot(plot_corr_T{:, vars{i}}, 'Positions', x_box, 'Colors', 'k',...
        'Symbol', '','Widths', 0.4, 'MedianStyle', 'line');

    % Loop over protocols and plot each swarm
    for j = 1:numProtocols
        protocol_j = protocols{j};
        idx = plot_corr_T.protocol == protocol_j;
        s = swarmchart(ones(sum(idx),1), plot_corr_T{idx, vars{i}}, ...
            100, repmat(colors(j,:), sum(idx),1), 'filled', ...
            'MarkerFaceAlpha', 0.2,'MarkerEdgeColor','k',...
            'MarkerEdgeAlpha',0.3);
        s.XJitterWidth = 0.3;
    end
    
    ylabel(varNames{i});
    xticks([]); % No x-ticks for individual subplots
    title(vars{i});
    xlim([0.5 1.5])
    ylim([0 (1.25*max(corr_T{:, vars{i}}))])

    % Add a legend
    lgd = legend({'Steps', 'Regular pulses', 'Variable pulses'}, 'Location', 'bestoutside');
    title(lgd, 'Protocol');

    % Set the figure size
    set(gcf, 'Position', [100, 100, 400, 600]);

    med = median(plot_corr_T{:, vars{i}});

    y_pos = max(plot_corr_T{:, vars{i}});  % Or choose a fixed y-position

    text(1, y_pos, sprintf('Median = %.5f', med), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 10, ...
    'Color', 'r');

    tit = varNames{i} + "rmoutliers_values";

    saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
    saveas(gcf, fullfile(saveDir,(tit)), 'png')

end

%% R2 vs. AP Latency

figure;

x_data = corr_T.AP_latency;
y_data = corr_T.R2;

% Remove NaN values
valid_idx = ~isnan(x_data) & ~isnan(y_data); % Keep only valid (non-NaN) indices
x_clean = x_data(valid_idx);
y_clean = y_data(valid_idx);

%%% Plot scatterplot
scatter(x_clean, y_clean, 100, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color, 'LineWidth', 0.25,...
    'MarkerFaceAlpha', 0.25)
hold on;
% Plot y = 0 line
yline(0, 'k', 'LineWidth', 1); % Black solid line
hold on;

%%% Plot regression line
% Fit a linear model (y = mx + b)
p = polyfit(x_clean, y_clean, 1); % 1st-degree polynomial (linear)
% Generate fitted y-values
y_fit = polyval(p, x_clean);
% Plot the regression line
h = plot(x_clean, y_fit, 'LineWidth', 1); % Regression line
set(h, 'Color', [0.5 0.5 0.5]); % Force dashed line
% Calculate correlation coefficient (R-value) and p-value
[R, P] = corrcoef(x_clean, y_clean);
R_value = R(1,2); % Extract correlation coefficient
p_value = P(1,2); % Extract p-value
% Annotate the graph with R-value and p-value
text(min(x_clean) + 0.10, min(y_clean) + 0.1, sprintf('R = %.2f, p = %.8f', R_value, p_value), ...
    'FontSize', 12);

% Customize the plot
xlabel('Latency (ms)');
ylabel('R2 of Decay Model Fit');
%title('Fidelity vs. Current Injection');
%legend('show');  % Show legend
grid off;  % Add grid
hold off;  % Release the hold
%xlim([0 max(x_clean)])
%ylim([0 1])

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "R2_latency";

%saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
%saveas(gcf, fullfile(saveDir,(tit)), 'png')

%% deltaA vs. AP Latency

figure;

x_data = corr_T.AP_latency(corr_T.R2 > 0.5);
y_data = corr_T.param2(corr_T.R2 > 0.5);

% Remove NaN values
valid_idx = ~isnan(x_data) & ~isnan(y_data); % Keep only valid (non-NaN) indices
x_clean = x_data(valid_idx);
y_clean = y_data(valid_idx);

%%% Plot scatterplot
scatter(x_clean, y_clean, 100, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color, 'LineWidth', 0.25,...
    'MarkerFaceAlpha', 0.25)
hold on;
% Plot y = 0 line
yline(0, 'k', 'LineWidth', 1); % Black solid line
hold on;

%%% Plot regression line
% Fit a linear model (y = mx + b)
p = polyfit(x_clean, y_clean, 1); % 1st-degree polynomial (linear)
% Generate fitted y-values
y_fit = polyval(p, x_clean);
% Plot the regression line
h = plot(x_clean, y_fit, 'LineWidth', 1); % Regression line
set(h, 'Color', [0.5 0.5 0.5]); % Force dashed line
% Calculate correlation coefficient (R-value) and p-value
[R, P] = corrcoef(x_clean, y_clean);
R_value = R(1,2); % Extract correlation coefficient
p_value = P(1,2); % Extract p-value
% Annotate the graph with R-value and p-value
text(min(x_clean) + 0.10, max(y_clean), sprintf('R = %.2f, p = %.8f', R_value, p_value), ...
    'FontSize', 12);

% Customize the plot
xlabel('Latency (ms)');
ylabel('deltaA value of Decay Model Fit');
%title('Fidelity vs. Current Injection');
%legend('show');  % Show legend
grid off;  % Add grid
hold off;  % Release the hold
%xlim([0 max(x_clean)])
%ylim([0 1])

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "deltaA_AP_latency";

%saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
%saveas(gcf, fullfile(saveDir,(tit)), 'png')

%% Tau vs. AP Latency

figure;

x_data = corr_T.AP_latency(corr_T.R2 > 0.5);
y_data = corr_T.param3(corr_T.R2 > 0.5);

% Remove NaN values
valid_idx = ~isnan(x_data) & ~isnan(y_data); % Keep only valid (non-NaN) indices
x_clean = x_data(valid_idx);
y_clean = y_data(valid_idx);

%%% Plot scatterplot
scatter(x_clean, y_clean, 100, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color, 'LineWidth', 0.25,...
    'MarkerFaceAlpha', 0.25)
hold on;
% Plot y = 0 line
yline(0, 'k', 'LineWidth', 1); % Black solid line
hold on;

%%% Plot regression line
% Fit a linear model (y = mx + b)
p = polyfit(x_clean, y_clean, 1); % 1st-degree polynomial (linear)
% Generate fitted y-values
y_fit = polyval(p, x_clean);
% Plot the regression line
h = plot(x_clean, y_fit, 'LineWidth', 1); % Regression line
set(h, 'Color', [0.5 0.5 0.5]); % Force dashed line
% Calculate correlation coefficient (R-value) and p-value
[R, P] = corrcoef(x_clean, y_clean);
R_value = R(1,2); % Extract correlation coefficient
p_value = P(1,2); % Extract p-value
% Annotate the graph with R-value and p-value
text(min(x_clean) + 0.10, max(y_clean), sprintf('R = %.2f, p = %.8f', R_value, p_value), ...
    'FontSize', 12);

% Customize the plot
xlabel('AP Latency (ms)');
ylabel('Tau value of Decay Model Fit');
%title('Fidelity vs. Current Injection');
%legend('show');  % Show legend
grid off;  % Add grid
hold off;  % Release the hold
%xlim([0 max(x_clean)])
%ylim([0 1])

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "tau_AP_latency";

saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
saveas(gcf, fullfile(saveDir,(tit)), 'png')

%% Tau vs. AP Latency (outliers removed)

figure;

x_data = corr_T.AP_latency(corr_T.R2 > 0.5 & corr_T.param3 < 199);
y_data = corr_T.param3(corr_T.R2 > 0.5 & corr_T.param3 < 199);

% Remove NaN values
valid_idx = ~isnan(x_data) & ~isnan(y_data); % Keep only valid (non-NaN) indices
x_clean = x_data(valid_idx);
y_clean = y_data(valid_idx);

%%% Plot scatterplot
scatter(x_clean, y_clean, 100, 'MarkerFaceColor', pv_color,...
    'MarkerEdgeColor', pv_color, 'LineWidth', 0.25,...
    'MarkerFaceAlpha', 0.25)
hold on;
% Plot y = 0 line
yline(0, 'k', 'LineWidth', 1); % Black solid line
hold on;

%%% Plot regression line
% Fit a linear model (y = mx + b)
p = polyfit(x_clean, y_clean, 1); % 1st-degree polynomial (linear)
% Generate fitted y-values
y_fit = polyval(p, x_clean);
% Plot the regression line
h = plot(x_clean, y_fit, 'LineWidth', 1); % Regression line
set(h, 'Color', [0.5 0.5 0.5]); % Force dashed line
% Calculate correlation coefficient (R-value) and p-value
[R, P] = corrcoef(x_clean, y_clean);
R_value = R(1,2); % Extract correlation coefficient
p_value = P(1,2); % Extract p-value
% Annotate the graph with R-value and p-value
text(min(x_clean) + 0.10, max(y_clean), sprintf('R = %.2f, p = %.8f', R_value, p_value), ...
    'FontSize', 12);

% Customize the plot
xlabel('AP Latency (ms)');
ylabel('Tau value of Decay Model Fit');
%title('Fidelity vs. Current Injection');
%legend('show');  % Show legend
grid off;  % Add grid
hold off;  % Release the hold
%xlim([0 max(x_clean)])
%ylim([0 1])

% Set the figure size
set(gcf, 'Position', [100, 100, 800, 600]);

tit = "tau_AP_latency_no_outliers";

saveas(gcf, fullfile(saveDir,(tit)), 'pdf')
saveas(gcf, fullfile(saveDir,(tit)), 'png')
%% Group by cell to visualize parameter coherence

num_colors = 100;
cmap = lines(num_colors); % Get a predefined colormap
randomized_cmap = cmap(randperm(num_colors), :); % Shuffle colors
colormap(randomized_cmap);

scatter(categorical(corr_T.ID), corr_T.param2, 100, categorical(corr_T.ID),...
    'filled', 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0.5)

xlabel('Cell Identity')
ylabel('deltaA Value from Model Fit')

figure;

colormap(randomized_cmap);

scatter(categorical(corr_T.ID), corr_T.param3, 100, categorical(corr_T.ID),...
    'filled', 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0.5)

xlabel('Cell Identity')
ylabel('Tau Value from Model Fit')