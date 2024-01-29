
clear CompMetrics 
% OPtFile = 'DK_2018_Optimize_revenue_Bins_loadN_noSorting_whole_year_hor';
% BaseFile = 'Baseline_DK_2018_spl';
% load (['D:\data\34_pettas\PhD\Optimization\Results\' OPtFile])
% Baseline = load (['D:\data\34_pettas\PhD\Optimization\Final\Baseline Results all\Spline\' BaseFile]);


%% Binning to values 
%%%%% try to bin the get the binned response here
% 1 bin data to wsp-V pairs
% 2 go through the TS tresponse and assign the damage and revenue value to each bin
% 3 sum the values in the bins bins 
% 4 
clear dat aa daty Binnedy
in.Output = Output;
 [Binnedy{1},dat{1}] = GetBinnedVals(Baseline);
 [Binnedy{2},dat{2}] = GetBinnedVals(in);


%% plotting
figure
Drev = 100*((Binnedy{2}.Rev./Binnedy{1}.Rev)-1);
Drev(isnan(Drev)) = 0;
surf(repmat(dat{1}.Prbins,length(dat{1}.Vbins),1),repmat(dat{1}.Vbins',1,length(dat{1}.Prbins)),Drev)
view(gca,[0 90]);
colormap(redblue(128))
colorbar
xlim([min(dat{1}.Prbins) max(dat{1}.Prbins)])
ylim([min(dat{1}.Vbins) max(dat{1}.Vbins)])
xlabel('Electricity price eur/MWh')
ylabel('Wind speed m/s')
title(['Revenue total ' num2str(100*(sum(sum(Binnedy{2}.Rev))/sum(sum(Binnedy{1}.Rev))-1),'%.2f') '%'])
clear Drev

figure
surf(repmat(dat{1}.Prbins,length(dat{1}.Vbins),1),repmat(dat{1}.Vbins',1,length(dat{1}.Prbins)),Binnedy{1, 2}.Pow  )
view(gca,[0 90]);
% colormap(redblue(128))
colorbar
xlim([min(dat{1}.Prbins) max(dat{1}.Prbins)])
ylim([min(dat{1}.Vbins) max(dat{1}.Vbins)])
xlabel('Electricity price eur/MWh')
ylabel('Wind speed m/s')
title(['Power per bin %'])


figure
Ddam = 100*((Binnedy{2}.DAM.BROop./Binnedy{1}.DAM.BROop)-1);
Ddam(isnan(Ddam)) = 0;
surf(repmat(dat{1}.Prbins,length(dat{1}.Vbins),1),repmat(dat{1}.Vbins',1,length(dat{1}.Prbins)),Ddam)
view(gca,[0 90]);
colormap(redblue(128))
colorbar
xlim([min(dat{1}.Prbins) max(dat{1}.Prbins)])
ylim([min(dat{1}.Vbins) max(dat{1}.Vbins)])
xlabel('Electricity price eur/MWh')
ylabel('Wind speed m/s')
title(['BROop total ' num2str(100*(sum(sum(Binnedy{2}.DAM.BROop))/sum(sum(Binnedy{1}.DAM.BROop))-1),'%.2f') '%'])
clear Ddam

figure
Ddam = 100*((Binnedy{2}.DAM.LSSTq./Binnedy{1}.DAM.LSSTq)-1);
Ddam(isnan(Ddam)) = 0;
surf(repmat(dat{1}.Prbins,length(dat{1}.Vbins),1),repmat(dat{1}.Vbins',1,length(dat{1}.Prbins)),Ddam)
view(gca,[0 90]);
colormap(redblue(128))
colorbar
xlim([min(dat{1}.Prbins) max(dat{1}.Prbins)])
ylim([min(dat{1}.Vbins) max(dat{1}.Vbins)])
xlabel('Electricity price eur/MWh')
ylabel('Wind speed m/s')
title(['LSSTq total ' num2str(100*(sum(sum(Binnedy{2}.DAM.LSSTq))/sum(sum(Binnedy{1}.DAM.LSSTq))-1),'%.2f') '%'])
clear Ddam

figure
Ddam = 100*((Binnedy{2}.DAM.TTMx./Binnedy{1}.DAM.TTMx)-1);
Ddam(isnan(Ddam)) = 0;
surf(repmat(dat{1}.Prbins,length(dat{1}.Vbins),1),repmat(dat{1}.Vbins',1,length(dat{1}.Prbins)),Ddam)
view(gca,[0 90]);
colormap(redblue(128))
colorbar
xlim([min(dat{1}.Prbins) max(dat{1}.Prbins)])
ylim([min(dat{1}.Vbins) max(dat{1}.Vbins)])
title(['TTMx total ' num2str(100*(sum(sum(Binnedy{2}.DAM.TTMx))/sum(sum(Binnedy{1}.DAM.TTMx))-1),'%.2f') '%'])
clear Ddam

figure
Ddam = 100*((Binnedy{2}.DAM.BRMz./Binnedy{1}.DAM.BRMz)-1);
Ddam(isnan(Ddam)) = 0;
surf(repmat(dat{1}.Prbins,length(dat{1}.Vbins),1),repmat(dat{1}.Vbins',1,length(dat{1}.Prbins)),Ddam)
view(gca,[0 90]);
colormap(redblue(128))
colorbar
xlim([min(dat{1}.Prbins) max(dat{1}.Prbins)])
ylim([min(dat{1}.Vbins) max(dat{1}.Vbins)])
xlabel('Electricity price eur/MWh')
ylabel('Wind speed m/s')
title(['BRMz total ' num2str(100*(sum(sum(Binnedy{2}.DAM.BRMz))/sum(sum(Binnedy{1}.DAM.BRMz))-1),'%.2f') '%'])
clear Ddam

figure
Ddam = 100*((Binnedy{2}.DAM.TBMy./Binnedy{1}.DAM.TBMy)-1);
Ddam(isnan(Ddam)) = 0;
surf(repmat(dat{1}.Prbins,length(dat{1}.Vbins),1),repmat(dat{1}.Vbins',1,length(dat{1}.Prbins)),Ddam)
view(gca,[0 90]);
colormap(redblue(128))
colorbar
xlim([min(dat{1}.Prbins) max(dat{1}.Prbins)])
ylim([min(dat{1}.Vbins) max(dat{1}.Vbins)])
xlabel('Electricity price eur/MWh')
ylabel('Wind speed m/s')
title(['TBMy total ' num2str(100*(sum(sum(Binnedy{2}.DAM.TBMy))/sum(sum(Binnedy{1}.DAM.TBMy))-1),'%.2f') '%'])
clear Ddam

% % figure
% % surf(repmat(dat.Prbins,length(dat.Vbins),1),repmat(dat.Vbins',1,length(dat.Prbins)),100*Binnedy.DAM.TBMy/sum(sum(Binnedy.DAM.TBMy)))
% % view(gca,[0 90]);
% % colorbar
% % xlim([min(dat.Prbins) max(dat.Prbins)])
% % ylim([min(dat.Vbins) max(dat.Vbins)])
% % title('TBMy')
% 


% figure;plot(Output.options.control.FunSettings.optimizer{1, 1}.Vbins ,10*Output.options.control.FunSettings.optimizer{1, 1}.Values(1,1:length(Output.options.control.FunSettings.optimizer{1, 1}.Vbins)),LineWidth=2)
% title('Power rating per wsp bin %')
% grid minor
% xlabel('Wind speed bin center m/s')
% ylabel('Power rating %')


% close all
figure
subplot(5,1,1)
plot(Output.Time{1,1},Baseline.Output.cum.Revenue/max(Baseline.Output.cum.Revenue),Output.Time{1,1},Output.cum.Revenue/max(Baseline.Output.cum.Revenue))
legend({'Baseline' 'Opt'})
title('Revenue')
grid minor


subplot(5,1,2)
plot (Output.Time{1,1},Baseline.Output.cum.DAM.BROop/max(Baseline.Output.cum.DAM.BROop),Output.Time{1,1},Output.cum.DAM.BROop/max(Baseline.Output.cum.DAM.BROop))
% legend({'Baseline' 'Opt'})
title('Damage BROop')
grid minor

subplot(5,1,3)
plot (Output.Time{1,1},Baseline.Output.cum.DAM.TBMy/max(Baseline.Output.cum.DAM.TBMy),Output.Time{1,1},Output.cum.DAM.TBMy/max(Baseline.Output.cum.DAM.TBMy))
% legend({'Baseline' 'Opt'})
title('Damage TBMy')
grid minor

subplot(5,1,4)
plot (Output.Time{1,1},Baseline.Output.inst.Prat*10,'o',Output.Time{1,1},Output.inst.Prat*10,'o')
% legend({'Baseline' 'Opt'})
title('Power level %')
legend({['Base shut % ',num2str(100*Baseline.Output.metrics.ShutDown.Perc) ]; ['Optcase shut % ',num2str(100*Output.metrics.ShutDown.Perc) ]})
grid minor

subplot(5,1,5)
plot (Output.Time{1,1},Baseline.Output.inst.IBC,'o',Output.Time{1,1},Output.inst.IBC,'o')
% legend({'Baseline' 'Opt'})
title('IBC activation')
legend({['Baseline % ',num2str(Baseline.Output.metrics.IBCactivation.Perc) ]; ['Optcase % ',num2str(Output.metrics.IBCactivation.Perc) ]})
grid minor


figure
loadChans = {'TBMx' 'TBMy' 'TBMz' 'BRMx' 'BRMy' 'BRMz' 'BROop' 'BRIp' 'TTMx' 'TTMy' 'TTMz' 'LSSMy' 'LSSMz' 'LSSTq' 'Revenue'};
for ii=1:(length(loadChans)-1)
    CompMetrics(ii) =100*( Output.cum.DAM.(loadChans{ii})(end)/Baseline.Output.cum.DAM.(loadChans{ii})(end) -1) ; %#ok<*EVLDOT>
end
CompMetrics(end+1) = 100*( Output.cum.Revenue(end)/Baseline.Output.cum.Revenue(end) -1) ;
plot(CompMetrics,'X','LineWidth',2,'MarkerSize',15)
grid on
set(gca,'YMinorGrid','on');
xlim([0 length(loadChans)+1])
% xlabel('Time')
ylabel('Delta %')
% legend(namesFil)
title('Deltas cumulative')
set(gca,'xtick',[1:length(loadChans)],'xticklabel',loadChans,'TickLabelInterpreter','none') %#ok<*NBRAK>



%% Functions 
function [Binnedy,dat] = GetBinnedVals(in)
Output = in.Output;
dat.Vedge = 3:1:25;
dat.Vbins = 3.5:1:24.5;
dat.dp= 5; %price step
dat.Indon = find(Output.V>=4 & Output.V<=24 & Output.Price>0 & Output.Price<1500 & ~isnan(Output.Price) );
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
          Binnedy.Pow(iV,iPr) = 0;
          Binnedy.DAM.BROop(iV,iPr) = 0;
          Binnedy.DAM.TBMy(iV,iPr) = 0;
          Binnedy.DAM.BRMz(iV,iPr) = 0;
          Binnedy.DAM.LSSTq(iV,iPr) = 0; 
          Binnedy.DAM.TTMx(iV,iPr) = 0; 
      else     
          curr.Ind1 = find(dat.binPr(daty.curIV)==iPr); % find the indices of the V bins (curIV is mapped to binV) 
          if isempty(curr.Ind1)
              Binnedy.Rev(iV,iPr) = 0;
              Binnedy.Pow(iV,iPr) = 0;
              Binnedy.DAM.BROop(iV,iPr) = 0;
              Binnedy.DAM.TBMy(iV,iPr) = 0;
              Binnedy.DAM.BRMz(iV,iPr) = 0;
              Binnedy.DAM.LSSTq(iV,iPr) = 0;
              Binnedy.DAM.TTMx(iV,iPr) = 0; 
          else
              curr.Ind2 = daty.curIV(curr.Ind1); % finding the indices corresponding to the operating only conditions
              aa  = dat.Indon(curr.Ind2);  % mapping the index back to the global index of the general TS. Brain fuck!
              Binnedy.Rev(iV,iPr) = sum(Output.inst.Revenue(aa,1));
              Binnedy.Pow(iV,iPr) = 10*mean(Output.inst.Prat(aa,1));
              Binnedy.DAM.BROop(iV,iPr) = sum(Output.inst.DAM.BROop(aa,1));
              Binnedy.DAM.TBMy(iV,iPr) = sum(Output.inst.DAM.TBMy(aa,1));
              Binnedy.DAM.BRMz(iV,iPr) = sum(Output.inst.DAM.BRMz(aa,1));
              Binnedy.DAM.LSSTq(iV,iPr) = sum(Output.inst.DAM.LSSTq(aa,1));
              Binnedy.DAM.TTMx(iV,iPr) = sum(Output.inst.DAM.TTMx(aa,1)); 
          end
      end
      clear aa curr
    end
    clear daty
end


end





