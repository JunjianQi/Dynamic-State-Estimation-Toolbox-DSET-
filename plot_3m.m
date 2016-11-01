% Copyright (C) 2016 Junjian Qi
% Jun. 10, 2015
% plot results for ekf and ukf
clear
close all
clc
load('KF_3m.mat');
posState = 1:12;
states = states(:,s_pos);

plotheight=20;
plotwidth=16;
subplotsx=2;
subplotsy=3;   
leftedge=1.15;
rightedge=0.2;
topedge=0.2;
bottomedge=2;
spacex=1.3;
spacey=0.6;
fontsize=5;    
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% setting the Matlab figure
f = figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

for i=1:subplotsx
    for ii=1:subplotsy
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        plot(tsequence,states(:,posState(i*3-(ii-1)))','-','Color',[0.5,0.5,0.5],'linewidth',3);
        hold on
        plot(tsequence,E_MM(posState(i*3-(ii-1)),:),'m--','linewidth',3);
        plot(tsequence,U_MM(posState(i*3-(ii-1)),:),'k--','linewidth',1.5);
        plot(tsequence,U_MM1(posState(i*3-(ii-1)),:),'-.','Color',[1,0.6,0],'linewidth',3);
        plot(tsequence,U_MM2(posState(i*3-(ii-1)),:),'r:','linewidth',1.5);
        plot(tsequence,U_MM3(posState(i*3-(ii-1)),:),':','Color',[0.5,0.5,1],'linewidth',3.5);
        plot(tsequence,U_MM_gps(posState(i*3-(ii-1)),:),'--','Color',[0,0.5,0],'linewidth',3);
        plot(tsequence,U_MM_sq(posState(i*3-(ii-1)),:),'-.','Color',[0.0,0.0,0.8],'linewidth',1.5);
        if ii==subplotsy
%             title(['Title (',num2str(i),')'])
        end
 
        if ii>1
            set(ax,'xticklabel',[])
        end
 
        if i>1
%             set(ax,'yticklabel',[])
        end
        if i==1
            ylabel(['Ylabel (',num2str(ii),')'])
        end
        if ii==1
           xlabel(['Ylabel (',num2str(i),')'])
        end
        if ii==1
            xlabel(['$\textrm{Time\,(second)}$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
        else
            xlabel('','FontName','Times New Roman','FontSize',28);
        end
        set(gca,'FontName','Times New Roman','fontsize',28);
        if i==1
            if ii==1
                ylabel(['$\delta_3$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 -1 25]);
            elseif ii==2
                ylabel(['$\delta_2$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 -1 25]);
            else
                ylabel(['$\delta_1$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 -1 25]);
            end
        elseif i==2
            if ii==1
                ylabel(['$\omega_3$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 376 381]);
            elseif ii==2
                ylabel(['$\omega_2$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 375 387]);
            else
                ylabel(['$\omega_1$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 372 382]);
            end
        elseif i==3
            if ii==1
                ylabel(['$e_{q3}^\prime$'],'interpreter','latex','FontName','Times New Roman','FontSize',18);
            elseif ii==2
                ylabel(['$e_{q2}^\prime$'],'interpreter','latex','FontName','Times New Roman','FontSize',18);
            else
                ylabel(['$e_{q1}^\prime$'],'interpreter','latex','FontName','Times New Roman','FontSize',18);
            end
        else
            if ii==1
                ylabel(['$e_{d3}^\prime$'],'interpreter','latex','FontName','Times New Roman','FontSize',18);
            elseif ii==2
                ylabel(['$e_{d2}^\prime$'],'interpreter','latex','FontName','Times New Roman','FontSize',18);
            else
                ylabel(['$e_{d1}^\prime$'],'interpreter','latex','FontName','Times New Roman','FontSize',18);
            end
            set(gca,'YTick',[1 1.5 2]);
        end
    end
end
h = legend('$\textrm{Real State}$','$\textrm{EKF Estimate}$','$\textrm{UKF-schol Estimate}$',...
    '$\textrm{UKF-}\kappa\,\textrm{Estimate}$','$\textrm{UKF-modified Estimate}$','$\textrm{UKF-}\Delta Q\,\textrm{Estimate}$',...
    '$\textrm{UKF-GPS Estimate}$','$\textrm{SR-UKF Estimate}$','Location','southeast');
set(h,'Interpreter','latex')
set(h,'Fontsize',22);
set(gcf,'Units','centimeters','Position',[10 1 34 24]);

%Saving eps with matlab and then marking pdf and eps and png with system commands
%filename=['test'];
%print(gcf, '-depsc2','-loose',[filename,'.eps']);
%system(['epstopdf ',filename,'.eps'])
%system(['convert -density 300 ',filename,'.eps ',filename,'.png'])

