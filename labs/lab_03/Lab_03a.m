clc
clear
close all

addpath('./Matlab_Functions')
addpath('./SAR_Data')

% ----------------------- LAB ON SPLINES - 1 ------------------------------

%% 1. SAR data 1D
% --- 1.1 Data import
% Import the .mat file containing the matlab variables related to the PS of
% the Citylife's skyscraper
load('CitylifePS.mat')

% --- 1.2 Create working directories
% They are required for geoSplinter software
mkdir('data_input')
mkdir('data_output')
mkdir('job')

% --- 1.3 Create the dataset input file
writematrix([day_rel', res2], ...    % residuals
    strcat('.', filesep, 'data_input', filesep, 'PS_citylife_res.txt'), 'Delimiter', 'space')

% --- 1.4 Make two holes in the dataset
% This is usefull when we have gap in our data (not this especific case)
h_in = 41;
h_fin = 60;
h_in2 = 200;
h_fin2 = 250;
displ_hole = res2;
day_rel_hole = day_rel;
day_rel_hole(h_in:h_fin) = [];
displ_hole(h_in:h_fin) = [];
day_rel_hole(h_in2:h_fin2) = [];
displ_hole(h_in2:h_fin2) = [];

% --- 1.5 Create the dataset input file with hole
writematrix([day_rel_hole', displ_hole], ...
    strcat('.', filesep, 'data_input', filesep, 'PS_citylife_res_hole.txt'), 'Delimiter', 'space')


%% 2. Examples with geoSplinter software
% ----------- 2.1 Reasonable number of splines -----------
suffix_1 = 'lin1';

% Job file creation
data_dim = 1;                                   % Dimension of dataset (1D/2D)
type_spl = 1;                                   % Type of splines (linear/cubic)
delta_d  = 1;                                   % Delta discretization (1delta/2delta)
file_inp = 'PS_citylife_res.txt';               % Input filename
file_out = strcat('PS_citylife_', suffix_1) ;   % Output filename
num_obs  = length(displ);                       % Number of observations
num_spl  = 24;                                  % Number of nodes (= splines)
t_in     = day_rel(1);                          % First abscissa (= time)
t_fin    = day_rel(end);                        % Last abscissa (= time)
lambda   = 0;                                   % Regularization parameter (λ)
num_sig  = 8;                                   % Number of significant digits

jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, ...
            num_spl, t_in, t_fin, lambda, num_sig, delta_d);

% Run geoSplinter_analysis with the job file
jobFile_execution(file_out)

% Visualize the results
[data_lin1, raster_lin1, par_lin1, sigmaPar_lin1] = ...
    geoSplinter(strcat('.', filesep, 'data_output', filesep, file_out), 'lin');


% ----------- 2.2 Low number of splines -----------
suffix_2 = 'lin2';

% Job file creation
data_dim = 1;                                   % Dimension of dataset (1D/2D)
type_spl = 1;                                   % Type of splines (linear/cubic)
delta_d  = 1;                                   % Delta discretization (1delta/2delta)
file_inp = 'PS_citylife_res.txt';               % Input filename
file_out = strcat('PS_citylife_', suffix_2) ;   % Output filename
num_obs  = length(displ);                       % Number of observations
num_spl  = 5;                                   % Number of nodes (= splines)
t_in     = day_rel(1);                          % First abscissa (= time)
t_fin    = day_rel(end);                        % Last abscissa (= time)
lambda   = 0;                                   % Regularization parameter (λ)
num_sig  = 8;                                   % Number of significant digits

jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, ...
            num_spl, t_in, t_fin, lambda, num_sig, delta_d);

% Run geoSplinter_analysis with the job file
jobFile_execution(file_out)

% Visualize the results
[data_lin2, raster_lin2, par_lin2, sigmaPar_lin2] = ...
    geoSplinter(strcat('.', filesep, 'data_output', filesep, file_out), 'lin');


% ----------- 2.3 Small regularization (λ = 10) -----------
suffix_3 = 'lin3';

% Job file creation
data_dim = 1;                                   % Dimension of dataset (1D/2D)
type_spl = 1;                                   % Type of splines (linear/cubic)
delta_d  = 1;                                   % Delta discretization (1delta/2delta)
file_inp = 'PS_citylife_res.txt';               % Input filename
file_out = strcat('PS_citylife_', suffix_3) ;   % Output filename
num_obs  = length(displ);                       % Number of observations
num_spl  = 24;                                  % Number of nodes (= splines)
t_in     = day_rel(1);                          % First abscissa (= time)
t_fin    = day_rel(end);                        % Last abscissa (= time)
lambda   = 10;                                  % Regularization parameter (λ)
num_sig  = 8;                                   % Number of significant digits

jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, ...
            num_spl, t_in, t_fin, lambda, num_sig, delta_d);

% Run geoSplinter_analysis with the job file
jobFile_execution(file_out)

% Visualize the results
[data_lin3, raster_lin3, par_lin3, sigmaPar_lin3] = ...
    geoSplinter(strcat('.', filesep, 'data_output', filesep, file_out), 'lin');


% ----------- 2.4 Huge regularization (λ = 1000) -----------
suffix_4 = 'lin4';

% Job file creation
data_dim = 1;                                   % Dimension of dataset (1D/2D)
type_spl = 1;                                   % Type of splines (linear/cubic)
delta_d  = 1;                                   % Delta discretization (1delta/2delta)
file_inp = 'PS_citylife_res.txt';               % Input filename
file_out = strcat('PS_citylife_', suffix_4) ;   % Output filename
num_obs  = length(displ);                       % Number of observations
num_spl  = 24;                                  % Number of nodes (= splines)
t_in     = day_rel(1);                          % First abscissa (= time)
t_fin    = day_rel(end);                        % Last abscissa (= time)
lambda   = 1000;                                % Regularization parameter (λ)
num_sig  = 8;                                   % Number of significant digits

jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, ...
            num_spl, t_in, t_fin, lambda, num_sig, delta_d);

% Run geoSplinter_analysis with the job file
jobFile_execution(file_out)

% Visualize the results
[data_lin4, raster_lin4, par_lin4, sigmaPar_lin4] = ...
    geoSplinter(strcat('.', filesep, 'data_output', filesep, file_out), 'lin');


% ----------- 2.5 Hole and regularization -----------
suffix_5 = 'lin5';

% Job file creation
data_dim = 1;                                   % Dimension of dataset (1D/2D)
type_spl = 1;                                   % Type of splines (linear/cubic)
delta_d  = 1;                                   % Delta discretization (1delta/2delta)
file_inp = 'PS_citylife_res_hole.txt';          % Input filename
file_out = strcat('PS_citylife_', suffix_5) ;   % Output filename
num_obs  = length(displ_hole);                  % Number of observations
num_spl  = 24;                                  % Number of nodes (= splines)
t_in     = day_rel(1);                          % First abscissa (= time)
t_fin    = day_rel(end);                        % Last abscissa (= time)
lambda   = 10;                                  % Regularization parameter (λ)
num_sig  = 8;                                   % Number of significant digits

jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, ...
            num_spl, t_in, t_fin, lambda, num_sig, delta_d);

% Run geoSplinter_analysis with the job file
jobFile_execution(file_out)

% Visualize the results
[data_lin5, raster_lin5, par_lin5, sigmaPar_lin5] = ...
    geoSplinter(strcat('.', filesep, 'data_output', filesep, file_out), 'lin');


% ----------- 2.5 Cubic splines with hole and reg. -----------
suffix_5 = 'cub1';

% Job file creation
data_dim = 1;                                   % Dimension of dataset (1D/2D)
type_spl = 2;                                   % Type of splines (linear/cubic)
file_inp = 'PS_citylife_res_hole.txt';          % Input filename
file_out = strcat('PS_citylife_', suffix_5) ;   % Output filename
num_obs  = length(displ_hole);                  % Number of observations
num_spl  = 12;                                  % Number of nodes (= splines)
t_in     = day_rel(1);                          % First abscissa (= time)
t_fin    = day_rel(end);                        % Last abscissa (= time)
lambda   = 10;                                  % Regularization parameter (λ)
num_sig  = 8;                                   % Number of significant digits

jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, ...
            num_spl, t_in, t_fin, lambda, num_sig);

% Run geoSplinter_analysis with the job file
jobFile_execution(file_out)

% Visualize the results
[data_cub1, raster_cub1, par_cub1, sigmaPar_cub1] = ...
    geoSplinter(strcat('.', filesep, 'data_output', filesep, file_out), 'cub');


%% 3. Splines interpolation of residuals from Fourier lab
% - we removed the trend using a polynomial
% - we estimated the harmonics to model the periodic behavior
% - now we interpolate the residual signal with the chosen number of 
%   splines and parameters

% Job file creation
data_dim = 1;                               % Dimension of dataset (1D/2D)
type_spl = 2;                               % Type of splines (linear/cubic)
file_inp = 'PS_citylife_res.txt';           % Input filename
file_out = strcat('PS_citylife_res_cub');   % Output filename
num_obs  = length(displ);                   % Number of observations
num_spl  = 12;                              % Number of nodes (= splines)
t_in     = day_rel(1);                      % First abscissa (= time)
t_fin    = day_rel(end);                    % Last abscissa (= time)
lambda   = 0;                               % Regularization parameter (λ)
num_sig  = 8;                               % Number of significant digits

jobFile_analysis(data_dim, type_spl, file_inp, file_out, num_obs, ...
            num_spl, t_in, t_fin, lambda, num_sig);

% Run geoSplinter_analysis with the job file
jobFile_execution(file_out)

% Visualize the results
[data_res, raster_res, par_res, sigmaPar_res] = ...
    geoSplinter(strcat('.', filesep, 'data_output', filesep, file_out), 'cub');


% Test on the model: χ²
% Is the splines modelling introducing any additional meaningful
% description of the signal? Let's perform a test on the residuals
% considering a constant model.
alpha = 0.05;

chi2_obs = std(res2)/1*(length(res2)-1);
chi2_lim = chi2inv(1-alpha, length(res2)-1);
fprintf('Test results, X^2_obs: %.4f (X^2_lim: %.4f)\n', abs(chi2_obs), chi2_lim);

% -> splines were not needed (statistically)!


%% 4. Results' interpretation
% Define the residuals of the final iteration
res3 = data_res(:,4);

% Compute the residuals' fft
res3_fft = fftshift(fft(res3));
res2_fft = fftshift(fft(res2));

% Plot
figure
plot(abs(res3_fft))
hold on
plot(abs(res2_fft))


% --- 4.1 Compare the residuals with the Fourier analysis
idx_cc = (size(displ,1)-1)/2+1;
res2_fft_red = res2_fft;

res2_fft_red(idx_cc-2) = 0;         % negative 1st peak
res2_fft_red(idx_cc-1) = 0;         % negative 2nd peak
res2_fft_red(idx_cc+1) = 0;         % positive 2nd peak
res2_fft_red(idx_cc+2) = 0;         % positive 1st peak

% Plot
figure
plot(t, res3)
hold on
plot(t, ifft(ifftshift(res2_fft_red)))


%% 5. Final signal model
% Sum the contributions from polynomial, Fourier and splines
mod_prev = displ - res2;
mod_fin = mod_prev + data_res(:,3);

% Plot
figure
plot(t, displ, 'k')
hold on
plot(t, mod_prev, 'b')
plot(t, mod_fin, 'r')