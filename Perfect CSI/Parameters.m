%% Description:
% This script defines the simulation parameters.

%% Simulation Parameters:
N=8;                    % Number of receive antennas at the BS
K=8;                   % Number of users connected to the BS
rho_dB=-10:1:20;        % Normalized transmit SNR [dB]
rho=10.^(rho_dB/10);    % Normalized transmit SNR
r=1;                  % Data rate [bpcu]
MC=1e5;                 % Number of Monte Carlo runs
rng(1);                 % Set the random number generator
