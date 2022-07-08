%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MRC-SIC Receiver
Pout=zeros(1,length(rho));
avg_SIC=zeros(1,length(rho));                               % Average number of SIC Operations
parfor x=1:length(rho)
    sigma2_e=1/(L*rho(x));                                  % Variance of the channel estimation errors
    E=sqrt(sigma2_e/2)*E_aux;                               % Matrix of channel estimation errors
    G_hat=G_true+E;                                         % Estimated channel matrix
    D=zeros(1,MC);                                          % Number of correctly decoded mMTC packets in each run    
    numSIC=zeros(1,MC);
    for j=1:MC
        G=G_true(:,:,j);                                    % True channel matrix
        A=G_hat(:,:,j);                                     % MRC linear detector matrix
        G_til=E(:,:,j);                                     % Matrix of channel estimation errors
        [~,ord]=sort(vecnorm(A,2,1).^2,'descend');          % Determine the SIC decoding ordering based on the SNR
        G=G(:,ord);                                         % Sort the matrix of true channel gains
        A=A(:,ord);                                         % Sort the matrix of estimated channel gains
        G_til=G_til(:,ord);                                 % Sort the matrix of channel estimation errors
        for k=1:K
            % Compute the interference from the MTC devices waiting to be
            % decoded:
            I1=0;
            for k1=k+1:K
                I1=I1+abs(A(:,k)'*G(:,k1))^2;                 
            end
            % Compute the residual interference from the MTC devices
            % previously decoded:
            I2=0;
            for k2=1:(k-1)
                I2=I2+abs(A(:,k)'*G_til(:,k2))^2;
            end
            % Compute the SINR of the k-th MTC device:
            sigma_k=(rho(x)*(abs(A(:,k)'*G(:,k))^2))/(rho(x)*I1+rho(x)*I2+(vecnorm(A(:,k)',2)^2));
            % If the decoding succeeds, update the number of correctly decoded MTC devices:
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
    Pout(x)=1-mean(D)/K;                   % Outage probability
    avg_SIC(x)=mean(numSIC);                    % Average Number of SIC Operations
end

%% Saving the results
Pout_MRC_SIC_iCSI=Pout;
avg_SIC_MRC_SIC_iCSI=avg_SIC;
save('Results_MRC_SIC_iCSI.mat','rho_dB','Pout_MRC_SIC_iCSI','avg_SIC_MRC_SIC_iCSI')

%% Plotting the results
fig1=figure(1);
    semilogy(rho_dB,Pout_MRC_SIC_iCSI,'LineWidth',1.5)    
    hold on
    grid on
    ylim([0 1])
    
fig2=figure(2);
    plot(rho_dB,avg_SIC_MRC_SIC_iCSI,'LineWidth',1.5)
    hold on
    grid on
    
%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end
