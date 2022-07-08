%% Reset
clearvars
% close all
clc

%% Loading the results
load('Results_MRC.mat')
load('Results_MRC_SIC.mat')
load('Results_ZF.mat')
load('Results_MMSE.mat')
load('Results_MMSE_SIC.mat')
load('Results_ZF_SIC.mat')

%% Plotting the results
fig1=figure(1);
    set(fig1,'Position',[300 300 500 420])  
    semilogy(rho_dB,Pr_Em_MRC,'--','LineWidth',1.5)
    hold
    semilogy(rho_dB,Pr_Em_MRC_SIC,'--','LineWidth',1.5)
    semilogy(rho_dB,Pr_Em_ZF,'--','LineWidth',1.5)
    semilogy(rho_dB,Pr_Em_MMSE,'--','LineWidth',1.5)
    semilogy(rho_dB,Pr_Em_ZF_SIC,'--','LineWidth',1.5)
    semilogy(rho_dB,Pr_Em_MMSE_SIC,'--','LineWidth',1.5)
    grid on
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    ylim([1e-4 1])
    xlabel('$\rho$ [dB]','Interpreter','latex','FontSize',12)
    ylabel('$\mathcal{P}_{\textrm{out}}$','Interpreter','latex','FontSize',12)
    leg=legend('MRC','MRC-SIC','ZF','MMSE','ZF-SIC','MMSE-SIC');
    set(leg,'Interpreter','latex','FontSize',12,'Location','northoutside','NumColumns',3)
    
% fig2=figure(2);
%     set(fig2,'Position',[300 300 500 420])  
%     plot(rho_dB,zeros(1,length(rho_dB)),'-*','LineWidth',1.5)
%     hold on
%     plot(rho_dB,zeros(1,length(rho_dB)),'-o','LineWidth',1.5)
%     plot(rho_dB,ones(1,length(rho_dB)),'-s','LineWidth',1.5)
%     plot(rho_dB,ones(1,length(rho_dB)),'-x','LineWidth',1.5)
%     plot(rho_dB,avg_Inv_ZF_SIC,'-d','LineWidth',1.5)
%     plot(rho_dB,avg_Inv_MMSE_SIC,'-+','LineWidth',1.5)
%     grid on
%     set(gca,'TickLabelInterpreter','latex','FontSize',12)
%     xlabel('$\rho$ [dB]','Interpreter','latex','FontSize',12)
%     ylabel('Average number of matrix inversions','Interpreter','latex','FontSize',12)
%     leg=legend('MRC','MRC-SIC','ZF','MMSE','ZF-SIC','MMSE-SIC');
%     set(leg,'Interpreter','latex','FontSize',12,'Location','northoutside','NumColumns',3)
%     
% fig3=figure(3);
%     set(fig3,'Position',[300 300 500 420])  
%     plot(rho_dB,zeros(1,length(rho_dB)),'-*','LineWidth',1.5)
%     hold on
%     plot(rho_dB,avg_SIC_MRC_SIC,'-o','LineWidth',1.5)
%     plot(rho_dB,zeros(1,length(rho_dB)),'-s','LineWidth',1.5)
%     plot(rho_dB,zeros(1,length(rho_dB)),'-x','LineWidth',1.5)
%     plot(rho_dB,avg_SIC_ZF_SIC,'-d','LineWidth',1.5)
%     plot(rho_dB,avg_SIC_MMSE_SIC,'-+','LineWidth',1.5)
%     grid on
%     set(gca,'TickLabelInterpreter','latex','FontSize',12)
%     xlabel('$\rho$ [dB]','Interpreter','latex','FontSize',12)
%     ylabel('Average number of SIC operations','Interpreter','latex','FontSize',12)
%     leg=legend('MRC','MRC-SIC','ZF','MMSE','ZF-SIC','MMSE-SIC');
%     set(leg,'Interpreter','latex','FontSize',12,'Location','northoutside','NumColumns',3)