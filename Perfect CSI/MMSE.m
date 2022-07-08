%% Reset
clearvars
% close all
clc

%% Load simulation parameters:
Parameters;
if K>N
   error('Error! For the MMSE receiver, K cannot be greater than N.'); 
end

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MMSE Receiver
Pout=zeros(1,length(rho));
for x=1:length(rho)
    D=zeros(1,MC);                                  % Number of correctly decoded devices in each run
    for j=1:MC
        G=G_MC(:,:,j);                              % Select the channel gains of the K devices
        temp=inv(eye(K)+rho(x)*G'*G);               % Temporary Matrix
        sigma_k=1./real(sum(diag(temp),2))-1;       % SINR of the mMTC devices being decoded
        D(j)=sum(log2(1+sigma_k)>=r);               % Update the number of correctly decoded devices
    end          
    Pout(x)=1-mean(D)/K;                            % Compute the outage probability
end

%% Saving the results
Pout_MMSE=Pout;
save('Results_MMSE.mat','rho_dB','Pout_MMSE')

%% Plotting the results
semilogy(rho_dB,Pout,'LineWidth',1.5)
grid on
ylim([0 1])

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end