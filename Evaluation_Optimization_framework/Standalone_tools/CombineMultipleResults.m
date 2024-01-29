clc 
clearvars
close all


% Script that takes reults and combines them to one to see the total
% effect. The files Areadded in the order defined


pathbase= 'D:\data\34_pettas\PhD\Optimization\Final\Baseline Results all\Spline\';
pathopt = 'D:\data\34_pettas\PhD\Optimization\Final\P_IBC_Shut_opt\';
OptFiles = {...
'DE_2016_REV';
'DE_2017_REV';
'DE_2018_REV';
 }';
BaseFiles = {
    'Baseline_DE_2016_spl'
    'Baseline_DE_2017_spl'
    'Baseline_DE_2018_spl'
    };

sensors= {
    'TBMx'
    'TBMy'
    'TBMz'
    'BRMx'
    'BRMy'
    'BRMz'
    'BROop'
    'BRIp'
    'TTMx'
    'TTMy'
    'TTMz'
    'LSSMy'
    'LSSMz'
    'LSSTq'
    };
for i = 1:length(OptFiles)
    Opt{i,1} = load([pathopt OptFiles{i}]);
    Base{i,1} = load([pathbase BaseFiles{i}]);
    Optcum.Energy(i,1) =  Opt{i}.Output.cum.Energy(end);
    Optcum.Rev(i,1) =  Opt{i}.Output.cum.Revenue(end) ;
    Optcum.IBCact(i,1) = Opt{i}.Output.metrics.IBCactivation.Perc ; 
    Optcum.Shutperc(i,1) =100*Opt{i}.Output.metrics.ShutDown.Perc  ;
    Basecum.Energy(i,1) =  Base{i}.Output.cum.Energy(end);
    Basecum.Rev(i,1) =  Base{i}.Output.cum.Revenue(end) ;
    Basecum.IBCact(i,1) = Base{i}.Output.metrics.IBCactivation.Perc ;
    Basecum.Shutperc(i,1) =100*Base{i}.Output.metrics.ShutDown.Perc  ;
    for iSen = 1:length(sensors)
        Optcum.(sensors{iSen})(i,1) = Opt{i}.Output.cum.DAM.(sensors{iSen})(end);
        Basecum.(sensors{iSen})(i,1) = Base{i}.Output.cum.DAM.(sensors{iSen})(end);
    end
end


for ii=1:(length(sensors))
    CompMetrics(ii) =100*(sum(Optcum.(sensors{ii}))/sum(Basecum.(sensors{ii}))-1);
end
CompMetrics(end+1) = 100*(sum(Optcum.Energy)/sum(Basecum.Energy)-1);
CompMetrics(end+1) = 100*(sum(Optcum.Rev)/sum(Basecum.Rev)-1);
CompMetrics(end+1) = mean(Optcum.IBCact);
% CompMetrics(end+1) = 100*(sum(Optcum.Shutperc*8760)/sum(Basecum.Shutperc*8760)-1);
plotleg = sensors;
plotleg{end+1} = 'Energy';
plotleg{end+1} = 'Revenue';
plotleg{end+1} = 'IBC on %';
% plotleg{end+1} = 'ShutDown';
figure
plot(CompMetrics,'X','LineWidth',2,'MarkerSize',15)
grid on
set(gca,'YMinorGrid','on');
xlim([0 length(plotleg)+1])
% xlabel('Time')
ylabel('Delta %')
% legend(namesFil)
title('Deltas to baseline cumulative combine DE 2016-2018')
set(gca,'xtick',[1:length(plotleg)],'xticklabel',plotleg,'TickLabelInterpreter','none') %#ok<*NBRAK>
legend(['shut Opt= ' num2str(mean(Optcum.Shutperc),'%2.1f') '% Baseline= ' num2str(mean(Basecum.Shutperc),'%2.1f') '%'],Location="southeast",FontSize=10)
aa= 1;