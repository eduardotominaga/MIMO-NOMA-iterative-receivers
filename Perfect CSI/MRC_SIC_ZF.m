%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MRC-SIC/ZF Receiver
Pout=zeros(1,length(rho));                              % Outage probability
avg_Inv=zeros(1,length(rho));                           % Average number of matrix inversions
avg_SIC=zeros(1,length(rho));                           % Average number of SIC operations
parfor x=1:length(rho)
    D=zeros(1,MC);                                      % Number of correctly decoded devices in each run
    num_Inv=zeros(1,MC);                                % Number of matrix inversions in each run
    num_SIC=zeros(1,MC);                                % Number of SIC operations in each run
    for j=1:MC
        G=G_MC(:,:,j);                                  % Select the channel gains of the K devices
        [~,ord]=sort(vecnorm(G,2,1).^2,'descend');
        G_ord=G(:,ord);                                 % Sort the devices in descending order of SNRs
        if K>N                      % If the number of devices is greater than the number of antennas...
            % MRC-SIC decoding :
            for k=1:(K-N)
                % Interference from the devices waiting to be decoded:
                I=0;
                for k1=k+1:K
                    I=I+abs(G_ord(:,k)'*G_ord(:,k1))^2;       
                end
                % SINR of the device being decoded:
                sigma_k=(rho(x)*vecnorm(G_ord(:,k),2,1)^4)/(rho(x)*I+vecnorm(G_ord(:,k),2,1)^2);
                if log2(1+sigma_k)>=r                   % If the decoding succeeds...
                    D(j)=k;                             % Update the number of correctly decoded mMTC packets
                    num_SIC(j)=k;                       % Update the number of SIC operations
                else
                    break;                              % Else, the SIC decoding procedure ends.
                end                               
            end
            % ZF decoding:
            if k==K-N                                          
                temp=inv(G(:,k+1:end)'*G(:,k+1:end));       % Temporary Matrix
                num_Inv(j)=num_Inv(j)+1;                    % Update the number of matrix inversions
                sigma_k=rho(x)./real(sum(diag(temp),2));    % SINR of the devices being decoded
                D(j)=D(j)+sum(log2(1+sigma_k)>=r);          % Update the number of correctly decoded mMTC packets
            end
        else                                                % Else, traditional ZF decoding is used
            temp=inv(G'*G);                                 % Temporary Matrix
            num_Inv(j)=num_Inv(j)+1;                        % Update the number of matrix inversions
            sigma_k=rho(x)./real(sum(diag(temp),2));        % SINR of the devices being decoded
            D(j)=sum(log2(1+sigma_k)>=r);                   % Update the number of correctly decoded devices
        end
    end          
    Pout(x)=1-mean(D)/K;                        % Compute the outage probability
    avg_Inv(x)=mean(num_Inv);                   % Compute the average number of matrix inversions
    avg_SIC(x)=mean(num_SIC);                   % Compute the average number of SIC operations
end

%% Saving the results
Pout_MRC_SIC_ZF=Pout;
avg_Inv_MRC_SIC_ZF=avg_Inv;
avg_SIC_MRC_SIC_ZF=avg_SIC;
save('Results_MRC_SIC_ZF.mat','rho_dB','Pout_MRC_SIC_ZF','avg_Inv_MRC_SIC_ZF','avg_SIC_MRC_SIC_ZF')

%% Plotting the results
fig1=figure(1);
    semilogy(rho_dB,Pout,'LineWidth',1.5)
    hold on
    grid on
    ylim([0 1])
    
fig2=figure(2);
    plot(rho_dB,avg_SIC,'LineWidth',1.5)
    hold on
    grid on
    
fig3=figure(3);
    plot(rho_dB,avg_Inv,'LineWidth',1.5)
    grid on

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="STO201028"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end