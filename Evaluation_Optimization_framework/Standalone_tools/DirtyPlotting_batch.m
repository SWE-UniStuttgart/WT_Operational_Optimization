
clearvars
close all
clc


pathbase= '..\Results\';
pathopt = '..\Results\';
OptFiles = {...
'DE_all_DayPrev_LOAD_constTSR.mat'
'DE_all_DayPrev_REV_constTSR.mat'
    }';

%%%%% Dirty patch to get all
% clear OptFiles
% % OptiflesAll=dir(fullfile(pathopt,'*DE_2016*mid*.mat'));
% OptiflesAll=dir(fullfile(pathopt,'*DK_2016*.mat'));
% for i =1:length(OptiflesAll)
%     OptFiles{i,1}= OptiflesAll(i).name;
% end


BaseFile = 'Baseline_DE_all_spl.mat';
for i=1:length(OptFiles)+1
    if i ==1
        baseLeg{1} = 'Baseline' ;
    else
        baseLeg{i} = OptFiles{i-1}; %#ok<*SAGROW> 
    end
end

% load files
for i =1:length(OptFiles)
    Opt{i,1} = load([pathopt OptFiles{i}]);
end
Baseline = load([pathbase BaseFile]);

% get binned values
for i =1:length(OptFiles)+1
    if i==1
    [ Binned{i}, dat{i}]  =  GetBinnedVals(Baseline);
    else
     [ Binned{i}, dat{i}] =  GetBinnedVals(Opt{i-1,1});
    end
end


loadChans = {'TBMx' 'TBMy' 'TBMz' 'BRMx' 'BRMy' 'BRMz' 'BROop' 'BRIp' 'TTMx' 'TTMy' 'TTMz' 'LSSMy' 'LSSMz' 'LSSTq' 'Energy' 'Revenue'};
plotChans ={'TBMx' 'TBMy' 'TBMz' 'BRMx' 'BRMy' 'BRMz' 'BROop' 'BRIp' 'TTMx' 'TTMy'  'LSSMy'  'LSSTq' };
% close all
figure
subplot(5,1,1)
plot(Opt{1}.Output.Time{1,1},Baseline.Output.cum.Revenue/max(Baseline.Output.cum.Revenue))
    hold on
for i =1:length(OptFiles)
    plot(Opt{i}.Output.Time{1,1},Opt{i}.Output.cum.Revenue/max(Baseline.Output.cum.Revenue))
    hold on
end
legend(baseLeg,Interpreter="none")
title('Revenue')
grid minor

subplot(5,1,2)
plot(Opt{i}.Output.Time{1,1},Baseline.Output.cum.DAM.BROop/max(Baseline.Output.cum.DAM.BROop))
hold on
for i =1:length(OptFiles)
    plot (Opt{i}.Output.Time{1,1},Opt{i}.Output.cum.DAM.BROop/max(Baseline.Output.cum.DAM.BROop))
    hold on
end
title('Damage BROop')
grid minor

subplot(5,1,3)
plot(Opt{1}.Output.Time{1,1},Baseline.Output.cum.DAM.TBMy/max(Baseline.Output.cum.DAM.TBMy))
hold on
for i =1:length(OptFiles)
    plot (Opt{i}.Output.Time{1,1},Opt{i}.Output.cum.DAM.TBMy/max(Baseline.Output.cum.DAM.TBMy))
    hold on
end
title('Damage TBMy')
grid minor

subplot(5,1,4)
plot(Opt{1}.Output.Time{1,1},Baseline.Output.inst.Prat*10,'o')
hold on
for i =1:length(OptFiles)
    plot (Opt{i}.Output.Time{1,1},Opt{i}.Output.inst.Prat*10,'o')
    hold on
end
title('Power level %')
legcell1{1,1} = ['Baseline % ',num2str(100*Baseline.Output.metrics.ShutDown.Perc,'%2.1f') ] ;
for i =1:length(OptFiles)
    legcell1{i+1,1} = ['Optcase shut % ',num2str(100*Opt{i}.Output.metrics.ShutDown.Perc,'%2.1f') ];
end
legend(legcell1)
grid minor

subplot(5,1,5)
plot(Opt{1}.Output.Time{1,1},Baseline.Output.inst.IBC,'o')
hold on
for i =1:length(OptFiles)
    plot (Opt{i}.Output.Time{1,1},Opt{i}.Output.inst.IBC,'o')
    %     legend({['Baseline % ',num2str(Baseline.Output.metrics.IBCactivation.Perc) ]; ['Optcase % ',num2str(Opt{i}.Output.metrics.IBCactivation.Perc) ]})
    hold on
end
title('IBC activation')
legcell2{1,1} = ['Baseline % ',num2str(Baseline.Output.metrics.IBCactivation.Perc,'%2.1f') ] ;
for i =1:length(OptFiles)
    legcell2{i+1,1} = ['Optcase IBC % ',num2str(Opt{i}.Output.metrics.IBCactivation.Perc,'%2.1f') ];
end
grid minor
legend(legcell2)


figure
for i =1:length(OptFiles)
    for ii=1:(length(loadChans)-2)
        CompMetrics(ii) =100*( Opt{i}.Output.cum.DAM.(loadChans{ii})(end)/Baseline.Output.cum.DAM.(loadChans{ii})(end) -1) ; %#ok<*EVLDOT>
    end
    
    CompMetrics(end+1) = 100*( Opt{i}.Output.cum.Energy(end)/Baseline.Output.cum.Energy(end) -1) ;
    CompMetrics(end+1) = 100*( Opt{i}.Output.cum.Revenue(end)/Baseline.Output.cum.Revenue(end) -1) ;
    plot(CompMetrics,'X','LineWidth',2,'MarkerSize',15)
    hold on
    clear CompMetrics
end
grid on
set(gca,'YMinorGrid','on');
% xlabel('Time')
ylabel('Delta %')
xlim([0 length(loadChans)+1])
legend(OptFiles,Interpreter="none",Location="southeast")
title('Deltas cumulative')
set(gca,'xtick',[1:length(loadChans)],'xticklabel',loadChans,'TickLabelInterpreter','none') %#ok<*NBRAK>


figure
for ii=1:length(plotChans)
        subplot(6,2,ii)
for i =1:length(OptFiles)
    plot (Opt{i}.Output.Time{1,1},Baseline.Output.cum.DAM.(plotChans{ii})/max(Baseline.Output.cum.DAM.(plotChans{ii})),Opt{i}.Output.Time{1,1},Opt{i}.Output.cum.DAM.(plotChans{ii})/max(Baseline.Output.cum.DAM.(plotChans{ii})))
    % legend({'Baseline' 'Opt'})
    hold on
end
title(plotChans{ii})
grid minor
end
legend(OptFiles,Interpreter="none")

% plot power per speed for constant prices
if isfield(Opt{1,1}.Output.options,'const_price') && Opt{1,1}.Output.options.const_price.activate  ==1
    figure
    for i =1:length(OptFiles)
        plot (Opt{i}.Output.options.control.FunSettings.optimizer{1, 1}.Vbins,Opt{i}.Output.options.control.FunSettings.optimizer{1,1}.Values(1:length(Opt{i}.Output.options.control.FunSettings.optimizer{1, 1}.Vbins)),LineWidth=2 )
        hold on
    end
    legend(OptFiles,Interpreter="none",Location="southeast")
    title('Power rating per wsp bin %')
    grid minor
    xlabel('Wind speed bin center m/s')
    ylabel('Power rating %')
end

%%% plot heat maps
% Revenue
for i =2:max(length(OptFiles),2)
    DRevbins = 100*(Binned{1,i}.Rev./Binned{1,1}.Rev-1)   ;
    DRevbins(isnan(DRevbins)) =0;
    figure
    surf(repmat(dat{i}.Prbins,length(dat{i}.Vbins),1),repmat(dat{i}.Vbins',1,length(dat{i}.Prbins)),DRevbins)
    view(gca,[0 90]);
    colormap(redblue(128))
    colorbar
    xlim([min(dat{i}.Prbins) max(dat{i}.Prbins)])
    ylim([min(dat{i}.Vbins) max(dat{i}.Vbins)])
    xlabel('Electricity price eur/MWh')
    ylabel('Wind speed m/s')
    title(['ΔRev case' num2str(i-1) ' total ' num2str(100*(sum(sum(Binned{i}.Rev))/sum(sum(Binned{1}.Rev))-1),'%.2f') '%'])
    clear DRevbins
end

% Lsstq
for i =2:max(length(OptFiles),2)
    DLSSTqbins = 100*(Binned{1,i}.DAM.LSSTq./Binned{1,1}.DAM.LSSTq-1)   ;
    DLSSTqbins(isnan(DLSSTqbins)) =0;
    figure
    surf(repmat(dat{i}.Prbins,length(dat{i}.Vbins),1),repmat(dat{i}.Vbins',1,length(dat{i}.Prbins)),DLSSTqbins)
    view(gca,[0 90]);
    colormap(redblue(128))
    colorbar
    xlim([min(dat{i}.Prbins) max(dat{i}.Prbins)])
    ylim([min(dat{i}.Vbins) max(dat{i}.Vbins)])
    xlabel('Electricity price eur/MWh')
    ylabel('Wind speed m/s')
    title(['ΔDAM LSSTQ case' num2str(i-1) ' total ' num2str(100*(sum(sum(Binned{i}.DAM.LSSTq))/sum(sum(Binned{1}.DAM.LSSTq))-1),'%.2f') '%'])
    clear DRevbins
end

% BROop
for i =2:max(length(OptFiles),2)
    DBROopbins = 100*(Binned{1,i}.DAM.BROop./Binned{1,1}.DAM.BROop-1)   ;
    DBROopbins(isnan(DBROopbins)) =0;
    figure
    surf(repmat(dat{i}.Prbins,length(dat{i}.Vbins),1),repmat(dat{i}.Vbins',1,length(dat{i}.Prbins)),DBROopbins)
    view(gca,[0 90]);
    colormap(redblue(128))
    colorbar
    xlim([min(dat{i}.Prbins) max(dat{i}.Prbins)])
    ylim([min(dat{i}.Vbins) max(dat{i}.Vbins)])
    xlabel('Electricity price eur/MWh')
    ylabel('Wind speed m/s')
    title(['ΔDAM BROop case' num2str(i-1) ' total ' num2str(100*(sum(sum(Binned{i}.DAM.BROop))/sum(sum(Binned{1}.DAM.BROop))-1),'%.2f') '%'])
    clear DRevbins
end


% TBMy
for i =2:max(length(OptFiles),2)
    DTBMybins = 100*(Binned{1,i}.DAM.TBMy./Binned{1,1}.DAM.TBMy-1)   ;
    DTBMybins(isnan(DTBMybins)) =0;
    figure
    surf(repmat(dat{i}.Prbins,length(dat{i}.Vbins),1),repmat(dat{i}.Vbins',1,length(dat{i}.Prbins)),DTBMybins)
    view(gca,[0 90]);
    colormap(redblue(128))
    colorbar
    xlim([min(dat{i}.Prbins) max(dat{i}.Prbins)])
    ylim([min(dat{i}.Vbins) max(dat{i}.Vbins)])
    xlabel('Electricity price eur/MWh')
    ylabel('Wind speed m/s')
    title(['ΔDAM TBMy case' num2str(i-1) ' total ' num2str(100*(sum(sum(Binned{i}.DAM.TBMy))/sum(sum(Binned{1}.DAM.TBMy))-1),'%.2f') '%'])
    clear DRevbins
end
%% Functions 
function [Binnedy,dat] = GetBinnedVals(in)
Output = in.Output;
dat.Vedge = 3:1:25;
dat.Vbins = 3.5:1:24.5;
dat.dp= 5; %price step
dat.Indon = find(Output.V>=4 & Output.V<=24 & Output.Price>0 & ~isnan(Output.Price) );
dat.V = Output.V(dat.Indon);
dat.Price = Output.Price(dat.Indon);
dat.Predge = 0:dat.dp:max(dat.Price)+dat.dp;
dat.Prbins= dat.dp/2:dat.dp:max(dat.Predge)-dat.dp/2;
[~,~,~,dat.binV,dat.binPr]=histcounts2(dat.V,dat.Price,dat.Vedge,dat.Predge);
for iV =1:length(dat.Vbins)
    daty.curIV =  find(dat.binV==iV);  
    for iPr = 1:length(dat.Prbins)     
      if isempty(daty.curIV)
          Binnedy.Rev(iV,iPr) = 0;
          Binnedy.DAM.BROop(iV,iPr) = 0;
          Binnedy.DAM.TBMy(iV,iPr) = 0;
          Binnedy.DAM.BRMz(iV,iPr) = 0;
          Binnedy.DAM.LSSTq(iV,iPr) = 0; 
      else     
          curr.Ind1 = find(dat.binPr(daty.curIV)==iPr); % find the indices of the V bins (curIV is mapped to binV) 
          if isempty(curr.Ind1)
              Binnedy.Rev(iV,iPr) = 0;
              Binnedy.DAM.BROop(iV,iPr) = 0;
              Binnedy.DAM.TBMy(iV,iPr) = 0;
              Binnedy.DAM.BRMz(iV,iPr) = 0;
              Binnedy.DAM.LSSTq(iV,iPr) = 0;
          else
              curr.Ind2 = daty.curIV(curr.Ind1); % finding the indices corresponding to the operating only conditions
              aa  = dat.Indon(curr.Ind2);  % mapping the index back to the global index of the general TS. Brain fuck!
              Binnedy.Rev(iV,iPr) = sum(Output.inst.Revenue(aa,1));
              Binnedy.DAM.BROop(iV,iPr) = sum(Output.inst.DAM.BROop(aa,1));
              Binnedy.DAM.TBMy(iV,iPr) = sum(Output.inst.DAM.TBMy(aa,1));
              Binnedy.DAM.BRMz(iV,iPr) = sum(Output.inst.DAM.BRMz(aa,1));
              Binnedy.DAM.LSSTq(iV,iPr) = sum(Output.inst.DAM.LSSTq(aa,1));
          end
      end
      clear aa curr
    end
    clear daty
end


end


