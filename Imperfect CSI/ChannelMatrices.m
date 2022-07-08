%% Description:
% This script generates the wireless channel matrices, including an
% auxiliary matrix for the generation of random channel estimation errors.

%% Wireless channel matrices:
G_true=(1/sqrt(2))*(randn(N,K,MC)+1i*randn(N,K,MC));    % Matriz of channel coefficients of all the mMTC Devices
E_aux=randn(N,K,MC)+1i*randn(N,K,MC);                   % Auxiliary matrix of channel estimation errors