%% Reset
clearvars
close all
clc

%% Loading the results

% Perfect CSI:
load('Results_MRC.mat')
load('Results_MRC_SIC.mat')
load('Results_ZF.mat')
load('Results_MMSE.mat')
load('Results_MMSE_PIC.mat')
load('Results_ZF_PIC.mat')

% Imperfect CSI:
load('Results_MRC_iCSI.mat')
load('Results_MRC_SIC_iCSI.mat')
load('Results_ZF_iCSI.mat')
load('Results_MMSE_iCSI.mat')
load('Results_MMSE_PIC_iCSI.mat')
load('Results_ZF_PIC_iCSI.mat')

%% Define the colors for the plots
blue = [57 106 177]./255;
red = [204 37 41]./255;
black = [83 81 84]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;
yellow = [240 180 50]./255;

%% Plotting the results

% Plot the outage probability versus transmit SNR:
fig1=figure(1);
    set(fig1,'Position',[300 300 600 450])  
    p1=semilogy(rho_dB,Pout_MRC,'-','LineWidth',2,'color',black);
    hold
    paux1=semilogy(rho_dB,Pout_MRC,'-','LineWidth',2,'color',black);
    p2=semilogy(rho_dB,Pout_MRC_iCSI,'--','LineWidth',2,'color',black);
    paux2=semilogy(rho_dB,Pout_MRC_iCSI,'--','LineWidth',2,'color',black);
    p3=semilogy(rho_dB,Pout_MRC_SIC,'-','LineWidth',2,'color',red);
    p4=semilogy(rho_dB,Pout_MRC_SIC_iCSI,'--','LineWidth',2,'color',red);
    p5=semilogy(rho_dB,Pout_ZF,'-','LineWidth',2,'color',blue);
    p6=semilogy(rho_dB,Pout_ZF_iCSI,'--','LineWidth',2,'color',blue);
    p7=semilogy(rho_dB,Pout_MMSE,'-','LineWidth',2,'color',green);
    p8=semilogy(rho_dB,Pout_MMSE_iCSI,'--','LineWidth',2,'color',green);
    p9=semilogy(rho_dB,Pout_ZF_PIC,'-','LineWidth',2,'color',yellow);
    p10=semilogy(rho_dB,Pout_ZF_PIC_iCSI,'--','LineWidth',2,'color',yellow);
    p11=semilogy(rho_dB,Pout_MMSE_PIC,'-','LineWidth',2,'color',purple);
    p12=semilogy(rho_dB,Pout_MMSE_PIC_iCSI,'--','LineWidth',2,'color',purple);
    grid on
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    ylim([1e-4 1])
    xlabel('$\rho$ [dB]','Interpreter','latex','FontSize',14)
    ylabel('$\mathcal{P}_{\textrm{out}}$','Interpreter','latex','FontSize',14)
%     leg1=legend([p1 p3 p5 p7 p9 p11 paux1 paux2],'MRC','MRC-SIC','ZF','MMSE','ZF-PIC','MMSE-PIC');
%     set(leg1,'Interpreter','latex','FontSize',14,'NumColumns',3,'Location','northoutside')
    a=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(a,[paux1 paux2],'Perfect CSI','Imperfect CSI');
    set(leg2,'Interpreter','latex','FontSize',14,'Location','southeast')
    saveas(fig1,'Figure_1A.png')
    saveas(fig1,'Figure_1A.eps','epsc')

% Plot the average number of matrix inversions versus transmit SNR:
fig2=figure(2);
    set(fig2,'Position',[300 300 600 450])  
    p1=plot(rho_dB,zeros(1,length(rho_dB)),'-','LineWidth',2,'color',black);       % MRC, perfect CSI    
    hold on
    paux1=plot(rho_dB,zeros(1,length(rho_dB)),'-','LineWidth',2,'color',black);    % MRC, perfect CSI
    p2=plot(rho_dB,zeros(1,length(rho_dB)),'--','LineWidth',2,'color',black);      % MRC, imperfect CSI
    paux2=plot(rho_dB,zeros(1,length(rho_dB)),'--','LineWidth',2,'color',black);   % MRC, imperfect CSI    
    p3=plot(rho_dB,zeros(1,length(rho_dB)),'-','LineWidth',2,'color',red);       % MRC-SIC, perfect CSI
    p4=plot(rho_dB,zeros(1,length(rho_dB)),'--','LineWidth',2,'color',red);      % MRC-SIC, imperfect CSI
    p5=plot(rho_dB,ones(1,length(rho_dB)),'-','LineWidth',2,'color',blue);        % ZF, perfect CSI
    p6=plot(rho_dB,ones(1,length(rho_dB)),'--','LineWidth',2,'color',blue);       % ZF, imperfect CSI
    p7=plot(rho_dB,ones(1,length(rho_dB)),'-','LineWidth',2,'color',green);        % MMSE, perfect CSI
    p8=plot(rho_dB,ones(1,length(rho_dB)),'--','LineWidth',2,'color',green);       % MMSE, imperfect CSI
    p9=plot(rho_dB,avg_Inv_ZF_PIC,'-','LineWidth',2,'color',yellow);
    p10=plot(rho_dB,avg_Inv_ZF_PIC_iCSI,'--','LineWidth',2,'color',yellow);
    p11=plot(rho_dB,avg_Inv_MMSE_PIC,'-','LineWidth',2,'color',purple);
    p12=plot(rho_dB,avg_Inv_MMSE_PIC_iCSI,'--','LineWidth',2,'color',purple);
    grid on
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    xlabel('$\rho$ [dB]','Interpreter','latex','FontSize',14)
    ylabel('Average number of matrix inversions','Interpreter','latex','FontSize',14)
%     leg=legend([p1 p3 p5 p7 p9 p11],'MRC','MRC-SIC','ZF','MMSE','ZF-PIC','MMSE-PIC');
%     set(leg,'Interpreter','latex','FontSize',14,'Location','northoutside','NumColumns',3)
    a=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(a,[paux1 paux2],'Perfect CSI','Imperfect CSI');
    set(leg2,'Interpreter','latex','FontSize',14,'Location','southeast')
    saveas(fig2,'Figure_1B.png')
    saveas(fig2,'Figure_1B.eps','epsc')

% Plot the average number of SIC operations versus transmit SNR:
fig3=figure(3);
    set(fig3,'Position',[300 300 600 550])  
    p1=plot(rho_dB,zeros(1,length(rho_dB)),'-','LineWidth',2,'color',black);       % MRC, perfect CSI    
    hold on
    paux1=plot(rho_dB,zeros(1,length(rho_dB)),'-','LineWidth',2,'color',black);    % MRC, perfect CSI
    p2=plot(rho_dB,zeros(1,length(rho_dB)),'--','LineWidth',2,'color',black);      % MRC, imperfect CSI
    paux2=plot(rho_dB,zeros(1,length(rho_dB)),'--','LineWidth',2,'color',black);   % MRC, imperfect CSI    
    p3=plot(rho_dB,avg_SIC_MRC_SIC,'-','LineWidth',2,'color',red);               % MRC-SIC, perfect CSI
    p4=plot(rho_dB,avg_SIC_MRC_SIC_iCSI,'--','LineWidth',2,'color',red);         % MRC-SIC, imperfect CSI
    p5=plot(rho_dB,zeros(1,length(rho_dB)),'-','LineWidth',2,'color',blue);       % ZF, perfect CSI
    p6=plot(rho_dB,zeros(1,length(rho_dB)),'--','LineWidth',2,'color',blue);      % ZF, imperfect CSI
    p7=plot(rho_dB,zeros(1,length(rho_dB)),'-','LineWidth',2,'color',green);       % MMSE, perfect CSI
    p8=plot(rho_dB,zeros(1,length(rho_dB)),'--','LineWidth',2,'color',green);      % MMSE, imperfect CSI
    p9=plot(rho_dB,avg_SIC_ZF_PIC,'-','LineWidth',2,'color',yellow);
    p10=plot(rho_dB,avg_SIC_ZF_PIC_iCSI,'--','LineWidth',2,'color',yellow);
    p11=plot(rho_dB,avg_SIC_MMSE_PIC,'-','LineWidth',2,'color',purple);
    p12=plot(rho_dB,avg_SIC_MMSE_PIC_iCSI,'--','LineWidth',2,'color',purple);
    grid on
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    xlabel('$\rho$ [dB]','Interpreter','latex','FontSize',14)
    ylabel('Average number of SIC operations','Interpreter','latex','FontSize',14)
    leg=legend([p1 p3 p5 p7 p9 p11],'MRC','MRC-SIC','ZF','MMSE','ZF-PIC','MMSE-PIC');
    set(leg,'Interpreter','latex','FontSize',14,'Location','southoutside','NumColumns',3)
    a=axes('position',get(gca,'position'),'visible','off');
    leg2=legend(a,[paux1 paux2],'Perfect CSI','Imperfect CSI');
    set(leg2,'Interpreter','latex','FontSize',14,'Location','East')
    saveas(fig3,'Figure_1C.png')
    saveas(fig3,'Figure_1C.eps','epsc')
    