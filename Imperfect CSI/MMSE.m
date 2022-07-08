%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;
if K>N
    error('Error! For the MMSE receiver, K cannot be greater than N.');
end

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MMSE Receiver
Pout=zeros(1,length(rho));                          % Outage probability
parfor x=1:length(rho)
    sigma2_e=1/(L*rho(x));                          % Variance of the channel estimation errors
    E=sqrt(sigma2_e/2)*E_aux;                       % Matrix of channel estimation errors
    G_hat=G_true+E;                                 % Estimated channel matrix
    D=zeros(1,MC);                                  % Number of correctly decoded MTC devices in each iteration
    for j=1:MC
        G=G_true(:,:,j);                                                    % Select the true channel matrix
        A=inv(G_hat(:,:,j)*G_hat(:,:,j)'+(1/rho(x))*eye(N))*G_hat(:,:,j);   % MMSE linear detector matrix
        for k=1:K
            % Compute the interference from other devices:
            I=0;
            for k1=1:K
                if k1~=k
                    I=I+abs(A(:,k)'*G(:,k1))^2;
                end                    
            end
            % SINR of the MTC device being decoded:
            sigma_k=(rho(x)*(abs(A(:,k)'*G(:,k))^2))/(rho(x)*I+(vecnorm(A(:,k)',2)^2));
            % If the device is correctly decoded, update the number of
            % correctly decoded devices:
            if log2(1+sigma_k)>=r
                D(j)=D(j)+1;
            end
        end
    end          
    Pout(x)=1-mean(D)/K;        % Compute the outage probability
end

%% Saving the results
Pout_MMSE_iCSI=Pout;
save('Results_MMSE_iCSI.mat','rho_dB','Pout_MMSE_iCSI')

%% Plotting the results
semilogy(rho_dB,Pout_MMSE_iCSI,'LineWidth',1.5)
grid on
ylim([0 1])

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end
