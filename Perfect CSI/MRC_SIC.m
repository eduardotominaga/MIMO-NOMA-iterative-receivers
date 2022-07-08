%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MRC-SIC Receiver
Pout=zeros(1,length(rho));                              % Outage probability
avg_SIC=zeros(1,length(rho));                           % Average number of SIC Operations
parfor x=1:length(rho)
    D=zeros(1,MC);                                      % Number of correctly decoded devices in each run
    numSIC=zeros(1,MC);                                 % Number of SIC operations in each run
    for j=1:MC
        G=G_MC(:,:,j);                                  % Select the channel gains of the K devices
        [~,ord]=sort(vecnorm(G,2,1).^2,'descend');
        G_ord=G(:,ord);                                 % Sort the devices in descending order of SNRs
        for k=1:K
            % Interference from the devices waiting to be decoded:
            I=0;
            for k1=k+1:K
                I=I+abs(G_ord(:,k)'*G_ord(:,k1))^2; 
            end
            % SINR of the device being decoded:
            sigma_k=(rho(x)*vecnorm(G_ord(:,k),2,1)^4)/(rho(x)*I+vecnorm(G_ord(:,k),2,1)^2);
            % If the decoding succeeds, update the number of correctly
            % decoded devices:
            if log2(1+sigma_k)>=r
                D(j)=k;
                if D(j)~=K          % If there are still devices waiting to be decoded...
                   numSIC(j)=k;     % The receiver performs a SIC operation.
                end
            else
                break;
            end
        end
    end          
    Pout(x)=1-mean(D)/K;                    % Compute the outage probability
    avg_SIC(x)=mean(numSIC);                % Compute the average Number of SIC Operations
end

%% Saving the results
Pout_MRC_SIC=Pout;
avg_SIC_MRC_SIC=avg_SIC;
save('Results_MRC_SIC.mat','rho_dB','Pout_MRC_SIC','avg_SIC_MRC_SIC')

%% Plotting the results
fig1=figure(1);
    semilogy(rho_dB,Pout_MRC_SIC,'LineWidth',1.5)
    hold on
    grid on
    ylim([0 1])   
    
fig2=figure(2);
    plot(rho_dB,avg_SIC_MRC_SIC,'LineWidth',1.5)
    hold on
    grid on

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end