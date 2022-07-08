%% Reset
clearvars
% close all
clc

%% Load simulation parameters:
Parameters;
if K>N
   error('Error! For the ZF receiver, K cannot be greater than N.'); 
end

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - ZF Receiver
Pout=zeros(1,length(rho));                              % Outage probability
avg_Inv=zeros(1,length(rho));                           % Average number of matrix inversions
parfor x=1:length(rho)
    D=zeros(1,MC);                                      % Number of correctly decoded devices in each run
    num_Inv=zeros(1,MC);                                % Number of matrix inversions in each run
    for j=1:MC
        G=G_MC(:,:,j);                                  % Select the channel gains of the K devices
        temp=inv(G'*G);                                 % Temporary Matrix
        num_Inv(j)=num_Inv(j)+1;                        % Update the number of matrix inversions
        sigma_k=rho(x)./real(sum(diag(temp),2));        % SINR of the K devices being decoded
        D(j)=sum(log2(1+sigma_k)>=r);                   % Update the number of correctly decoded mMTC packets
    end          
    Pout(x)=1-mean(D)/K;                                % Compute the outage probability
    avg_Inv(x)=mean(num_Inv);                           % Compute the average number of matrix inversions
end

%% Saving the results
Pout_ZF=Pout;
avg_Inv_ZF=avg_Inv;
save('Results_ZF.mat','rho_dB','Pout_ZF','avg_Inv_ZF')

%% Plotting the results
fig1=figure(1);
    semilogy(rho_dB,Pout,'LineWidth',1.5)
    grid on
    ylim([0 1])   

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end
