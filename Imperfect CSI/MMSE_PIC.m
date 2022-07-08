%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - MMSE-PIC Receiver
Pout=zeros(1,length(rho));                                  % Outage probability
avg_Inv=zeros(1,length(rho));                               % Average number of matrix inversions
avg_SIC=zeros(1,length(rho));                               % Average number of SIC operations
parfor x=1:length(rho)
    sigma2_e=1/(L*rho(x));                                  % Variance of the channel estimation errors
    E=sqrt(sigma2_e/2)*E_aux;                               % Matrix of channel estimation errors
    G_hat=G_true+E;                                         % Estimated channel matrix
    D=zeros(1,MC);                                          % Number of decoded users in each run
    num_Inv=zeros(1,MC);                                    % Number of matrix inversions in each run
    num_SIC=zeros(1,MC);                                    % Number of SIC operations in each run
    for j=1:MC
        G=G_true(:,:,j);                                    % Select the matrix of true channel gains
        G_h=G_hat(:,:,j);                                   % Select the matrix of estimated channel gains
        G_til=E(:,:,j);                                     % Matrix of channel estimation errors
        [~,ord]=sort(vecnorm(G,2,1).^2,'descend');          % Determine the SIC decoding ordering based on the SNR
        G=G(:,ord);                                         % Sort the matrix of true channel gains
        G_h=G_h(:,ord);                                     % Sort the matrix of estimated channel gains
        G_til=G_til(:,ord);                                 % Sort the matrix of channel estimation errors
        buffer_CSI_errors=[];                               % Initialize a buffer matrix of CSI errors
        flag=0;        
        while D(j)<K
            % Initialization necessary for the first iteration of the SIC
            % decoding procedure:
            if flag==0
                new_G=G;
                new_G_h=G_h;
            end
            if size(new_G,2)>N              % If there are more remaining devices than receive antennas...
                G_mmse=new_G(:,1:N);
                G_h_mmse=new_G_h(:,1:N);
            else                            % Else, if there are more receive antennas than remaining devices...
                G_mmse=new_G;
                G_h_mmse=new_G_h;
            end                
            A=inv(G_h_mmse*G_h_mmse'+(1/rho(x))*eye(N))*G_h_mmse;   % MMSE linear detector matrix
            num_Inv(j)=num_Inv(j)+1;                                % Update the number of matrix inversions
            sigma=zeros(1,size(G_mmse,2));                          % Initialize a vector for the SINR values
            for k=1:size(G_mmse,2)
                % Compute the interference from the MTC devices waiting to be
                % decoded:
                I1=0;
                for k1=1:size(new_G,2)
                    if k1~=k
                        I1=I1+abs(A(:,k)'*new_G(:,k1))^2; 
                    end                        
                end
                % Compute the residual interference from the devices
                % previously decoded:
                I2=0;
                for k2=1:size(buffer_CSI_errors,2)
                    I2=I2+abs(A(:,k)'*buffer_CSI_errors(:,k2))^2;                        
                end
                % Compute the SINR of the devices being simultaneously
                % decoded:
                sigma(k)=(rho(x)*(abs(A(:,k)'*G_mmse(:,k))^2))/(rho(x)*I1+rho(x)*I2+(vecnorm(A(:,k)',2)^2));
            end
            indexDec=log2(1+sigma)>=r;                      % Index of decoded devices
            indexNotDec=~indexDec;                          % Index of not decoded devices 
            numDec=sum(indexDec);                           % Compute the number of correctly decoded devices
            if numDec==0
                break                                       % If no device was correctly decoded, the procedure ends
            else
                D(j)=D(j)+numDec;                               % Update the number of correctly decoded devices
                if D(j)~=K                                      % If there still are devices waiting to be decoded...
                    num_SIC(j)=num_SIC(j)+numDec;               % Update the number of SIC operations
                end 
                G_failed=new_G(:,indexNotDec);                  % True channel vectors of the not decoded devices
                G_h_failed=new_G_h(:,indexNotDec);              % Estimated channel vectors of the not decoded devices
                buffer_CSI_errors=[buffer_CSI_errors G_mmse(:,indexDec)-G_h_mmse(:,indexDec)];    % Update the buffer of CSI errors
                new_G=[G_failed new_G(:,N+1:end)];              % Update the matrix of true channel gains for the next iteration
                new_G_h=[G_h_failed new_G_h(:,N+1:end)];        % Update the matrix of estimated channel gains for the next iteration                                
                flag=1;     
            end                                                                                                                                    
        end    
    end
    Pout(x)=1-mean(D)/K;                % Compute the outage probability
    avg_Inv(x)=mean(num_Inv);           % Compute the average number of matrix inversions
%     avg_SIC(x)=mean(D);                 % Compute the average number of SIC operations
    avg_SIC(x)=mean(num_SIC);           % Compute the average number of SIC operations
end

%% Saving the Results
Pout_MMSE_PIC_iCSI=Pout;
avg_Inv_MMSE_PIC_iCSI=avg_Inv;
avg_SIC_MMSE_PIC_iCSI=avg_SIC;
save('Results_MMSE_PIC_iCSI.mat','rho_dB','Pout_MMSE_PIC_iCSI','avg_Inv_MMSE_PIC_iCSI','avg_SIC_MMSE_PIC_iCSI')

%% Plotting the results
fig1=figure(1);
    semilogy(rho_dB,Pout_MMSE_PIC_iCSI,'LineWidth',1.5)
    hold on
    grid on
    ylim([0 1])
    
fig2=figure(2);
    plot(rho_dB,avg_SIC_MMSE_PIC_iCSI,'LineWidth',1.5)
    hold on
    grid on

fig3=figure(3);
    plot(rho_dB,avg_Inv_MMSE_PIC_iCSI,'LineWidth',1.5)
    grid on

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end
