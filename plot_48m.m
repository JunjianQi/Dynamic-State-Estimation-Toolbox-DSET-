% Copyright (C) 2016 Junjian Qi
% plot results for ukf-gps and sr-ukf

clear
close all
clc

flag = 1; % 1 for positive; 2 for srukf
if flag==1
    load('KF_48m.mat');
    mac_tra_idx = para{7};
    states = states(:,s_pos);
end

plotheight=20;
plotwidth=16;
subplotsx=2;
subplotsy=2;   
leftedge=1.69;
rightedge=0.4;   
topedge=1;
bottomedge=2.5;
spacex=1.6;
spacey=0.8;
fontsize=5;
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% set Matlab figure
f = figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
 
% choose the method
X = U_MM;     % UKF-schol
% X = U_MM_gps; % UKF-GPS
% X = U_MM_sq;  % SR-UKF
% X = U_MM1;    % UKF-kappa
% X = U_MM2;    % UKF-modified
% X = U_MM3;    % UKF-DeltaQ

for i=1:subplotsx
    for ii=1:subplotsy
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        hold on
        if ii==subplotsy
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
                ylabel(['$e_q^{\prime}$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                hold on
                for k=1:27
                    plot(tsequence,states(:,2*n_mac+k)-X(2*n_mac+k,:)','LineWidth',1.5); %,'k','LineWidth',1.5);
                end            
                axis([-0.1 10.1 -1 1.4]); % UKF
                axis([-0.1 10.1 -0.15 0.07]);
            elseif ii==2
                ylabel(['$\delta$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 -7000 7000]); % UKF
                axis([-0.1 10.1 -1 1]);
                hold on
                for k=1:48
                    plot(tsequence,states(:,k)-X(k,:)','LineWidth',1.5); %'k','LineWidth',1.5);
                end
            else
            end
%             axis([0 5 0.44 0.52]);
%             set(gca,'YTick',[0.15 0.16]);
        elseif i==2
            if ii==1
                ylabel(['$e_d^{\prime}$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                hold on
                for k=1:27
                    plot(tsequence,states(:,2*n_mac+27+k)-X(2*n_mac+27+k,:)','LineWidth',1.5); %'k','LineWidth',1.5);
                end   
                axis([-0.1 10.1 -2 2.5]); % UKF
                axis([-0.1 10.1 -0.45 0.2]);
            elseif ii==2
                ylabel(['$\omega$'],'interpreter','latex','FontName','Times New Roman','FontSize',28);
                axis([-0.1 10.1 -3500 3500]); % UKF
                axis([-0.1 10.1 -8 8]);
                hold on
                for k=1:48
                    plot(tsequence,states(:,k+48)-X(k+48,:)','LineWidth',1.5); % 'k','LineWidth',1.5);
                end
            else
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

% set(gcf,'Units','centimeters','Position',[10 -2 22 15]);
set(gcf,'Units','centimeters','Position',[10 -1 36 24]);

%Saving eps with matlab and then marking pdf and eps and png with system commands
%filename=['test'];
%print(gcf, '-depsc2','-loose',[filename,'.eps']);
%system(['epstopdf ',filename,'.eps'])
%system(['convert -density 300 ',filename,'.eps ',filename,'.png'])

