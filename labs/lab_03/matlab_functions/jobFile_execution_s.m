function jobFile_execution_s(filename)

% Function to automatically execute the job file of geoSplinter on both 
% Windows and macOS.
% 
% INPUT VARIABLES:
% Job file name: filename
%
%
% Roberto Monti
% Politecnico di Milano
% Last update: November 2024


% Job file execution
if ismac
    job_execution = strcat('./geoSplinter_synthesis_macOS < ./job/', filename, '.job');
elseif ispc
    job_execution = strcat('.\geoSplinter_synthesis.exe < .\job\', filename, '.job');
end
    
system(job_execution)
system('exit')