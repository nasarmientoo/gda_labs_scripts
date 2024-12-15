clc
clear
close all

addpath('./SAR_Data')

% --------------------- LAB ON FOURIER ANALYSIS ---------------------------

%% 1. SAR data
% ! Same as Lab_01.m !
% Import the csv file with the displacement time series for the PS
% representing a skyscraper in Milan citylife
data = readmatrix('PS_citylife.csv');

% Extract time, coordinates and displacement
displ = data(2, 4:end)';
day_rel = data(1, 4:end);
coords = data(2, 1:2);

% Create a datetime array for the dates
t0 = datetime(2018, 01, 06, 00, 00, 00);
t = t0 + day_rel;

% Plot the time series
figure
plot(t, displ)


%% 1. Find a polynomial trend and remove it from the data
% Convert dates into relative days
t_d = (day_rel - mean(day_rel))';

% LS solution
Ap = [t_d, t_d.^2, t_d.^3];
xp = inv(Ap' * Ap) * Ap' * displ;
yp = Ap * xp;

vp = displ - yp;
s02 = vp' * vp / (size(Ap,1) - size(Ap,2));
Cxx = s02 * inv(Ap' * Ap);

% Plot the interpolation
figure
plot(t, displ)
hold on
plot(t, yp)

% Check the significance of each parameter of the polynomial - t-test
alpha = 0.05;
t_obs = xp ./ sqrt(diag(Cxx));
t_lim = tinv(1-alpha/2, size(Ap,1) - size(Ap,2));

% --- Iteration 1 - recompute the LS solution
Ap = [t_d, t_d.^2];
xp = inv(Ap' * Ap) * Ap' * displ;
yp = Ap * xp;

vp = displ - yp;
s02 = vp' * vp / (size(Ap,1) - size(Ap,2));
Cxx = s02 * inv(Ap' * Ap);

% t-test
t_obs = xp ./ sqrt(diag(Cxx));
t_lim = tinv(1-alpha/2, size(Ap,1) - size(Ap,2));

% --- Iteration 2 - recompute the LS solution
Ap = t_d.^2;
xp = inv(Ap' * Ap) * Ap' * displ;
yp = Ap * xp;

vp = displ - yp;
s02 = vp' * vp / (size(Ap,1) - size(Ap,2));
Cxx = s02 * inv(Ap' * Ap);
% t-test
t_obs = xp ./ sqrt(diag(Cxx));
t_lim = tinv(1-alpha/2, size(Ap,1) - size(Ap,2));

% Plot the final interpolation
figure
plot(t,displ)
hold on
plot(t, yp)

% Remove the estimated polynomial
displ1 = displ - yp;


%% 2. Matlab functions for Fourier analysis/synthesis
% Understand how the Matlab built-in functions for the Fourier analysis and 
% synthesis work.

% Create a column vector with 5 values
x = ([1 2 1 3 1.5]');

% Compute the DFT (Discrete Fourier Transform)
% The fft function estimates the Fourier coefficients with the constant
% term as first value of the estimation vector (X1). So, the first half of
% the spectrum will be made by the first (n-1)/2 elements.
X0 = fft(x);
f0 = [0; 1; 2; -2; -1];

    % Print
    fprintf('%10s %20s\n', 'Frequency', 'Fcoeff');
    for i = 1:length(f0)
        fprintf('%10d %20s\n', f0(i), num2str(X0(i)));
    end

% For plotting and interpretation reasons it is easier and better to put
% the constant term (corresponding to zero frequency) to the middle of the
% spectrum so that it is symmetric. To reach this goal, the last (n-1)/2 
% elements are mirrored to becomes the first (n-1)/2 elements.
% Example:
% X1(1) -> X((n-1)/2)     = X(3)
% X1(2) -> X((n-1)/2 + 1) = X(4)
% X1(3) -> X((n-1)/2 + 2) = X(5)
% X1(4) -> X((n-1)/2 - 2) = X(1)
% X1(5) -> X((n-1)/2 - 1) = X(2)
% We can achieve this using the fftshift function as follows.
X = fftshift(fft(x));
f1 = [-2;-1;0;1;2];
f11 = (-(5-1)/2 : 1 : (5-1)/2)';
    
    % Print
    fprintf('%10s %20s\n', 'Frequency', 'Fcoeff');
    for i = 1:length(f1)
        fprintf('%10d %20s\n', f1(i), num2str(X(i)));
    end

% The matlab function ifft can be used to perform the Fourier synthesis,
% which means retriveing the signal from the values of the coefficients of
% the harmonics. Of course, if the order of the parameters has been altered
% through the fftshift function, than the same re-arrangement of the
% elements inside the vector has to be made.
x0 = ifft(X0);
disp(x0);

x1 = ifft(ifftshift(X));
disp(x1);

% Plot the Fourier amplitude spectrum
figure
stem(f1, abs(X));
xlabel('Frequency')
ylabel('Amplitude')

% Find the relationship between lecture-given and matlab formulas
a0 = mean(x);
a1 = mean(x.*cos(2*pi*(0:1:4)'/5*1));
a2 = mean(x.*cos(2*pi*(0:1:4)'/5*2));
b1 = mean(x.*sin(2*pi*(0:1:4)'/5*1));
b2 = mean(x.*sin(2*pi*(0:1:4)'/5*2));
X_man = [a0 0; a1 b1; a2 b2];


%% 2a. Even number of data
% How to compute the frequencies in case the number of data is even.
z = [x; 4.1];
Z = fftshift(fft(z));

fZ = (-6/2 : 1 : (6/2 -1))';

    % Print
    fprintf('%10s %20s\n', 'Frequency', 'Fcoeff');
    for i = 1:length(fZ)
        fprintf('%10d %20s\n', fZ(i), num2str(Z(i)));
    end


%% 3. White noise
% Understand what is the meaning of the white noise and how its Fourier 
% spectrum looks like.

% Creation of a random vector with num elements
num = 1001;
wn = randn(1, num);

% Compute the frequency of each element
% The frequencies are symmetric with respect to 0 (origin of the axis)
% and both positive and negative side will be half of the total number of
% data (this is valid only if the data is odd).
freq_wn = (-(num-1)/2 : 1 : (num-1)/2)';

% Fourier analysis: compute the Fourier coefficients of the white noise
% vector
cF_wn = fftshift(fft(wn));

% Plot the Fourier amplitude spectrum of the white noise
figure
plot(freq_wn, abs(cF_wn));
xlabel('Frequency')
ylabel('Amplitude')


%% 4. Computation of the Fourier Spectrum of the SAR data
% The dominant period is the inverse of the frequency: f = 1/T.

% Compute the frequencies for the spectrum
nObs = length(displ1);
freq_PS = 1/(6 * nObs) * (-(nObs-1)/2 : 1 : (nObs-1)/2)';

% Compute the Fourier coefficients
cF_1 = fftshift(fft(displ1));

% Plot
figure
plot(freq_PS, abs(cF_1));
xlabel('Frequency')
ylabel('Amplitude')


%% 5. Evaluation of the peaks in the spectrum
% --- Find the index of the peak frequency (= fundamental frequency)
[f_max, idx_max] = max(abs(cF_1));
idx_cc = (size(cF_1,1)-1)/2+1;   % index of the constant term (central frequency (=0))

% Plot for manual research of the index
figure
plot(abs(cF_1));

% Creation of an empty vector (all frequencies amplitudes are zero) apart
% from the peaks. Remember that the spectrum is symmetric.
cF_peak = zeros(size(cF_1));
cF_peak(146) = cF_1(146);         % negative peak
cF_peak(156) = cF_1(156);         % positive peak
cF_peak(idx_cc) = cF_1(idx_cc);   % constant

% Plot the "reduced" spectrum
figure
stem(freq_PS, abs(cF_peak))
xlabel('Frequency')
ylabel('Amplitude')


%% 6. Fourier syntesis
% Estimation of the signal in time corresponding to the extracted peaks.
% This means that a single sinusoidal harmonic will decribe the
% seasonality of our signal.

% Perform the synthesis (IFFT)
displ1_harm = ifft(ifftshift(cF_peak));

% Plot the synthesis and the raw time series
figure
plot(t, displ1)
hold on
plot(t, displ1_harm, 'LineWidth', 2)
xlabel('Time')
ylabel('Displacement [mm]')
legend('Raw time series', 'Fourier synthesis T_f')


%% 7. Error estimation
% Error estimate based on a significance test on s02. X^2 test to compare
% the variance with the one of the data.

% Residuals
res = displ1 - displ1_harm;

% Compute the precision of the estimated Fourier coefficients (const. + peak)
s02 = res' * res / (nObs - 3);

% Variance of the data
s02_th = 1;

% X^2 test
chi0_2_obs = s02 / s02_th * (nObs - 3);
chi0_2_lim = chi2inv(1-alpha/2, (nObs - 3));
fprintf('Peak harmonic, X^2_obs: %.4f (X^2_lim: %.4f)\n', abs(chi0_2_obs), chi0_2_lim);

%% 8. Two harmonics
% As the previous X^2 test failed, let's try adding the second dominant
% harmonic.

% Find the index of the first two peak frequencies
[~, idx_max2] = maxk(abs(cF_1), 4);

% Empty vector apart from the peaks.
cF_peak2 = zeros(size(cF_1));
cF_peak2(146) = cF_1(146);         % negative 1st peak
cF_peak2(156) = cF_1(156);         % positive 1st peak
cF_peak2(149) = cF_1(149);         % negative 2nd peak
cF_peak2(153) = cF_1(153);         % positive 2nd peak
cF_peak2(idx_cc) = cF_1(idx_cc);   % constant

% Plot the "reduced" spectrum
figure
stem(freq_PS, abs(cF_peak2))
xlabel('Frequency')
ylabel('Amplitude')

% Synthesis
displ1_harm2 = ifft(ifftshift(cF_peak2));

% Plot the synthesis and the raw time series
figure
plot(t, displ1)
hold on
plot(t, displ1_harm, 'LineWidth', 2)
plot(t, displ1_harm2, 'LineWidth', 2)
xlabel('Time')
ylabel('Displacement [mm]')
legend('Raw time series', 'Fourier synthesis T_f', 'Fourier synthesis T_p_2')

% New residuals
res2 = displ1 - displ1_harm2;
s02 = res2' * res2 / (nObs - 5);

% X^2 test
chi0_2_obs = s02 / s02_th * (nObs - 5);
chi0_2_lim = chi2inv(1-alpha/2, (nObs - 5));
fprintf('2 Peak harmonics, X^2_obs: %.4f (X^2_lim: %.4f)\n', abs(chi0_2_obs), chi0_2_lim);

% Plot of residuals and its spectrum
figure
subplot(1,2,1)
plot(t, res2)
subplot(1,2,2)
plot(freq_PS, abs(cF_1)-abs(cF_peak2))


% --- Analytical expression of s02
s02_f = (sum(displ1.^2) - 1/nObs*cF_peak2(idx_cc)^2 - 1/nObs*2*abs(cF_peak2(150))^2 ...
    - 1/nObs*2*abs(cF_peak2(146))^2) / (nObs - 5);
