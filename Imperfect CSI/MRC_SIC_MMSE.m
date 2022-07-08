%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MRC-SIC/MMSE Receiver
Pout=zeros(1,length(rho));                                  % Outage probability
avg_Inv=zeros(1,length(rho));                               % Average number of matrix inversions
avg_SIC=zeros(1,length(rho));                               % Average number of SIC operations
parfor x=1:length(rho)
    sigma2_e=1/(L*rho(x));                                  % Variance of the channel estimation errors
    E=sqrt(sigma2_e/2)*E_aux;                               % Matrix of channel estimation errors
    G_hat=G_true+E;                                         % Estimated channel matrix
    D=zeros(1,MC);                                          % Number of correctly decoded MTC devices in each run
    num_Inv=zeros(1,MC);                                    % Number of matrix inversions in each run
    num_SIC=zeros(1,MC);                                    % Number of SIC operations
    for j=1:MC
        G=G_true(:,:,j);                                    % Select the channel gains of the MTC devices
        A_mrc=G_hat(:,:,j);                                 % MRC linear detector matrix
        G_til=E(:,:,j);                                     % Matrix of channel estimation errors
        [~,ord]=sort(vecnorm(A_mrc,2,1).^2,'descend');      % Define the SIC decoding ordering according to the SNRs
        G=G(:,ord);                                         % Sort the matrix of true channel gains
        A_mrc=A_mrc(:,ord);                                 % Sort the matrix of estimated channel gains
        G_til=G_til(:,ord);                                 % Sort the matrix of channel estimation errors
        if K>N                                              % If there are more devices than receive antenna elements...
            for k=1:(K-N)                                   % Start of the MRC-SIC decoding procedure
                % Compute the interference from the MTC devices waiting to
                % be decoded:
                I1=0;
                for k1=k+1:K
                    I1=I1+abs(A_mrc(:,k)'*G(:,k1))^2;
                end
                % Compute the residual interference from the MTC devices
                % waiting to be decoded:
                I2=0;
                for k2=1:(k-1)
                    I2=I2+abs(A_mrc(:,k)'*G_til(:,k2))^2;
                end                
                % Compute the SINR of the k-th MTC device:
                sigma_k=(rho(x)*(abs(A_mrc(:,k)'*G(:,k))^2))/(rho(x)*I1+rho(x)*I2+(vecnorm(A_mrc(:,k)',2)^2));
                if log2(1+sigma_k)>=r   % If the k-th MTC device is correctly decoded...
                    D(j)=k;                                 % Update the number of correctly decoded mMTC packets
                    num_SIC(j)=k;                           % Update the number of SIC operations
                else
                    break;              % Else, the SIC decoding procedure ends.
                end                               
            end                                            
            if k==K-N            % If the number of remaining users equals the number of antennas, use MMSE decoding:                
                G=G(:,k+1:end);                             % Select the channel gains from the remaining devices 
                A_mmse=inv(A_mrc(:,k+1:end)*A_mrc(:,k+1:end)'+(1/rho(x))*eye(N))*A_mrc(:,k+1:end);  % MMSE linear detector matrix
                num_Inv(j)=num_Inv(j)+1;                    % Update the number of matrix inversions
                % Attempt to simultaneously decode of each one of the remaining users using MMSE:
                for k3=1:N
                    % Interference from other active MTC devices:
                    I1=0;
                    for k1=1:N
                        if k1~=k3
                            I1=I1+abs(A_mmse(:,k3)'*G(:,k1))^2;       
                        end                    
                    end
                    % Compute the residual interference from the devices
                    % previously decoded using MRC-SIC:
                    I2=0;
                    for k2=1:(k-1)
                        I2=I2+abs(A_mmse(:,k3)'*G_til(:,k2))^2;
                    end  
                    % SINR of the MTC device being decoded:
                    sigma_k=(rho(x)*(abs(A_mmse(:,k3)'*G(:,k3))^2))/(rho(x)*I1+rho(x)*I2+(vecnorm(A_mmse(:,k3)',2)^2));
                    if log2(1+sigma_k)>=r       % If the MTC device is correctly decoded...
                        D(j)=D(j)+1;            % Update the number of correctly decoded MTC devices
                    end
                end
            end
        else        % If there are more receive antennas than devices, traditional MMSE decoding is used:
            G=G_true(:,:,j);                                                    % Select the true channel matrix
            A=inv(G_hat(:,:,j)*G_hat(:,:,j)'+(1/rho(x))*eye(N))*G_hat(:,:,j);   % MMSE linear detector matrix
            num_Inv(j)=num_Inv(j)+1;                                            % Update the number of matrix inversions
            % Attempt to simultaneously decode the K devices using MMSE:
            for k=1:K
                % Compute the interference from other MTC devices:
                I1=0;
                for k1=1:K
                    if k1~=k
                        I1=I1+abs(A(:,k)'*G(:,k1))^2;
                    end                    
                end
                % SINR of the MTC device being decoded:
                sigma_k=(rho(x)*(abs(A(:,k)'*G(:,k))^2))/(rho(x)*I1+(vecnorm(A(:,k)',2)^2));
                if log2(1+sigma_k)>=r       % If the device is correctly decoded...
                    D(j)=D(j)+1;            % Update the number of correctly decoded MTC packets
                end
            end
        end
    end          
    Pout(x)=1-mean(D)/K;                        % Compute the outage probability
    avg_Inv(x)=mean(num_Inv);                   % Compute the average number of matrix inversions
    avg_SIC(x)=mean(num_SIC);                   % Compute the average number of SIC operations
end

%% Saving the results
Pout_MRC_SIC_MMSE_iCSI=Pout;
avg_Inv_MRC_SIC_MMSE_iCSI=avg_Inv;
avg_SIC_MRC_SIC_MMSE_iCSI=avg_SIC;
save('Results_MRC_SIC_MMSE_iCSI.mat','rho_dB','Pout_MRC_SIC_MMSE_iCSI','avg_Inv_MRC_SIC_MMSE_iCSI','avg_SIC_MRC_SIC_MMSE_iCSI')

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
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end

