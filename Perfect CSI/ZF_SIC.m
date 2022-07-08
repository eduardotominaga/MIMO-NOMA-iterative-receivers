%% Reset
clearvars
close all
clc

%% Load simulation parameters:
Parameters;

%% Load wireless channel matrices:
ChannelMatrices;

%% Monte Carlo Simulation - ZF-SIC Receiver
Pout=zeros(1,length(rho));                                  % Outage probability
avg_Inv=zeros(1,length(rho));                               % Average number of matrix inversions
avg_SIC=zeros(1,length(rho));                               % Average number of SIC operations
parfor x=1:length(rho)
    disp(rho(x))
    D=zeros(1,MC);                                          % Number of decoded devices in each run
    num_Inv=zeros(1,MC);                                    % Number of matrix inversions in each run
    num_SIC=zeros(1,MC);                                    % Number of SIC operations in each run
    for j=1:MC
        G=G_MC(:,:,j);                                      % Select the channel vectors of the K devices
        [~,ord]=sort(vecnorm(G,2,1).^2,'descend');          % Define the SIC decoding ordering
        G=G(:,ord);                                         % Sort the devices in descending order of SNRs
        flag=0;
        while D(j)<K                    % While there are still devices waiting to be decoded...
            if flag==0                  % On the first iteration of the SIC decoding procedure...
                new_G=G(:,:);           % Select the channle vectors of the K devices
            end
            if size(new_G,2)>N          % If the number of remaining devices is greater than the number of antennas...
                G_zf=new_G(:,1:N);      % Select a NxN matrix
            else                        % Else...
                G_zf=new_G;             % Utilize the channel vectors of all the remaining devices
            end  
            A=G_zf*inv(G_zf'*G_zf);             % Compute the ZF linear detector matrix
            num_Inv(j)=num_Inv(j)+1;            % Update the number of matrix inversions
            
            % Interference from the devices waiting to be decoded:
            I=0;
            for k1=2:size(new_G,2)
                I=I+abs(A(:,1)'*new_G(:,k1))^2;                       
            end
            % SINR of the device being decoded:
            sigma=(rho(x)*(abs(A(:,1)'*G_zf(:,1))^2))/(rho(x)*I+(vecnorm(A(:,1)',2)^2));
            % Test if the device is correctly decoded:         
            if log2(1+sigma)>=r
                D(j)=D(j)+1;                % Update the number of correctly decoded devices
                if D(j)~=K                      % If there still are devices waiting to be decoded...
                    num_SIC(j)=num_SIC(j)+1;    % Update the number of SIC operations
                end 
                new_G=new_G(:,2:end);       % Update the channel matrix for the next iteration
                flag=1;
            else
                break;
            end        
        end
                
    end
    Pout(x)=1-mean(D)/K;                % Compute the outage probability
    avg_Inv(x)=mean(num_Inv);           % Compute the average number of matrix inversions
    avg_SIC(x)=mean(num_SIC);           % Compute the average number of SIC operations
end

%% Saving the Results
Pout_ZF_SIC=Pout;
avg_Inv_ZF_SIC=avg_Inv;
avg_SIC_ZF_SIC=avg_SIC;
save('Results_ZF_SIC.mat','rho_dB','Pout_ZF_SIC','avg_Inv_ZF_SIC','avg_SIC_ZF_SIC')

%% Plotting the results
fig1=figure(1);
    semilogy(rho_dB,Pout_ZF_SIC,'LineWidth',1.5)
    hold on
    grid on
    ylim([0 1])
    
fig2=figure(2);
    plot(rho_dB,avg_SIC_ZF_SIC,'LineWidth',1.5)
    hold on
    grid on

fig3=figure(3);
    plot(rho_dB,avg_Inv_ZF_SIC,'LineWidth',1.5)
    grid on

%% This part of the code terminates all the Matlab processes is the script run on a server:
if getenv('COMPUTERNAME')~="OY2106111"  % If this is not my personal computer...    
    exit;                               % Terminate all the Matlab processes
end
