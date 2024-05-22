clc;
clear all;
close all;
load('Capacity_K_1.mat')
load('Capacity_K_2.mat')
load('Capacity_K_5.mat')
load('Capacity_K_10.mat')
load('Capacity_max.mat')
load('NMSE_K_1.mat')
load('NMSE_K_2.mat')
load('NMSE_K_5.mat')
load('NMSE_K_10.mat')
figure
plot(NMSE_K_1,'r-.','linewidth', 3)
hold on
plot(NMSE_K_2,'g:','linewidth', 3)
hold on
plot(NMSE_K_5,'b--','linewidth', 3)
hold on
plot(NMSE_K_10,'m-','linewidth', 3)
legend(' K = 1',' K = 2',' K = 5',' K = 10','location','best')
xlabel('Number of transmit metasurface layers, L');
ylabel('Fitting NMSE, \Delta');
set(gca,'fontsize',14)


figure
plot(Capacity_K_1,'r-.','linewidth',3)
hold on
plot(Capacity_K_2,'g:','linewidth',3)
hold on
plot(Capacity_K_5,'b--','linewidth',3)
hold on
plot(Capacity_K_10,'m-','linewidth',3)
hold on
plot(Capacity_max,'k-.','linewidth',3)
legend(' K = 1',' K = 2',' K = 5',' K = 10',' Optimal','location','best')
xlabel('Number of transmit metasurface layers, L');
ylabel('Channel capacity, C [bps/Hz]');
set(gca,'fontsize',14)














