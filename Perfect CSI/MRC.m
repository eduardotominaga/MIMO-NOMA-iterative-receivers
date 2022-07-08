%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MRC Receiver
Pout=zeros(1,length(rho));          % Outage probability
parfor x=1:length(rho)
    D=zeros(1,MC);                  % Number of correctly decoded devices in each run
    for j=1:MC
        G=G_MC(:,:,j);              % Select the channel gains of the K devices
        for k=1:K
            % Compute the interference from other active devices:
            I=0;
            for k1=1:K
                if k1~=k
                    I=I+abs(G(:,k)'*G(:,k1))^2;
                end                    
            end
            % SINR of the k-th device:
            sigma_k=(rho(x)*vecnorm(G(:,k),2,1)^4)/(rho(x)*I+vecnorm(G(:,k),2,1)^2);                                       
            % If the decoding succeeds, update the number of correctly
            % decoded devices:
            if log2(1+sigma_k)>=r
                D(j)=D(j)+1;                             
            end
        end
    end          
    Pout(x)=1-mean(D)/K;            % Compute the outage probability
end

%% Saving the results
Pout_MRC=Pout;
save('Results_MRC.mat','rho_dB','Pout_MRC')

%% Plotting the results
semilogy(rho_dB,Pout,'LineWidth',1.5)
hold on
grid on
ylim([0 1])

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end
