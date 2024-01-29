clc
clearvars
close all

%Script to get the stat files and create the turbine repsonse 3d response
%over WSP, TI and power output of controllers
%
% TODOs
% -Add mean value except DEL for all load channels

%% INPUT

folderCnt = 'D:\data\34_pettas\PhD\Simulations\All_constTSR\'; % general folder path with all the results of a controller

%Define the dimensions to be read
TIvec = 2:2:24;
Pvec= 5:0.5:13;
Vvec= 4:1:24;
Poff = 0; % power in kW to consider the turbine off. Mainly used to avoid the weird results in low wsp/power

% Define Controller name
CntrNameInit = 'constTSR';
savename =[ CntrNameInit '_surrogate.mat' ];

%%
kk =0;
for iTI = 1:length(TIvec)
    kk=kk+1;
    disp(['Case ' num2str(kk)])
    for iV = 1:length(Vvec)

        % if IPC and below rated go the the folder with the non-IPC results
        if Vvec(iV)<11 && strcmp(CntrNameInit(end-3:end),'_IPC')
            CntrName = CntrNameInit(1:end-4);
        else
            CntrName = CntrNameInit;
        end
        for iP = 1:length(Pvec)
            %             tic
            % if power boosting all are the same look at the constTSR folder
            if Pvec(iP)>= 10
                folderN =  'D:\data\34_pettas\PhD\Simulations\All_constTSR\';
                if Vvec(iV)>11 && strcmp(CntrNameInit(end-3:end),'_IPC')
                    CntrName = 'constTSR_IPC';
                else
                    CntrName = 'constTSR';
                end
            else
                folderN =  folderCnt;
            end

            curSubFold1 = ['All_turb_' CntrName '_SD1_TI' num2str(TIvec(iTI),'%02.0f') '\Stats2\'];
            curfileN1 = [CntrName '_P' num2str(Pvec(iP)) '_WSP' num2str(Vvec(iV),'%02.0f')  '_TI' num2str(TIvec(iTI),'%02.0f') '_SD1' '_results_stats'  ];
            curfileN1 = regexprep(curfileN1, '\.', 'd');
            curfileN1 = [folderN curSubFold1 curfileN1]; %#ok<*AGROW>

            curSubFold2 = ['All_turb_' CntrName '_SD2_TI' num2str(TIvec(iTI),'%02.0f') '\Stats2\'];
            curfileN2 = [CntrName '_P' num2str(Pvec(iP)) '_WSP' num2str(Vvec(iV),'%02.0f')  '_TI' num2str(TIvec(iTI),'%02.0f') '_SD2' '_results_stats'  ];
            curfileN2 = regexprep(curfileN2, '\.', 'd');
            curfileN2 = [folderN curSubFold2 curfileN2]; %#ok<*AGROW>

            curSubFold3 = ['All_turb_' CntrName '_SD3_TI' num2str(TIvec(iTI),'%02.0f') '\Stats2\'];
            curfileN3 = [CntrName '_P' num2str(Pvec(iP)) '_WSP' num2str(Vvec(iV),'%02.0f')  '_TI' num2str(TIvec(iTI),'%02.0f') '_SD3' '_results_stats'  ];
            curfileN3 = regexprep(curfileN3, '\.', 'd');
            curfileN3 = [folderN curSubFold3 curfileN3]; %#ok<*AGROW>

            DataIN{1} = load(curfileN1);
            DataIN{2} = load(curfileN2);
            DataIN{3} = load(curfileN3);

            DataCnt.GenSpeedSTD.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'GenSpeed','Std'},DataIN{2}.Stats{'GenSpeed','Std'},DataIN{3}.Stats{'GenSpeed','Std'}]);
            DataCnt.GenSpeedSTD.std(iV,iTI,iP)  = std([DataIN{1}.Stats{'GenSpeed','Std'},DataIN{2}.Stats{'GenSpeed','Std'},DataIN{3}.Stats{'GenSpeed','Std'}]);

            DataCnt.Power.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'GenPwr','Mean'},DataIN{2}.Stats{'GenPwr','Mean'},DataIN{3}.Stats{'GenPwr','Mean'}]);
            DataCnt.Power.std(iV,iTI,iP)  = std([DataIN{1}.Stats{'GenPwr','Mean'},DataIN{2}.Stats{'GenPwr','Mean'},DataIN{3}.Stats{'GenPwr','Mean'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.Power.mean(iV,iTI,iP) = 0;
                DataCnt.Power.std(iV,iTI,iP) = 0;
            end

            DataCnt.Energy.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'GenPwr','Integr'},DataIN{2}.Stats{'GenPwr','Integr'},DataIN{3}.Stats{'GenPwr','Integr'}]);
            DataCnt.Energy.std(iV,iTI,iP) = std([DataIN{1}.Stats{'GenPwr','Integr'},DataIN{2}.Stats{'GenPwr','Integr'},DataIN{3}.Stats{'GenPwr','Integr'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.Energy.mean(iV,iTI,iP) = 0;
                DataCnt.Energy.std(iV,iTI,iP) = 0;
            end


            DataCnt.GenTqSTD.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'GenTq','Std'},DataIN{2}.Stats{'GenTq','Std'},DataIN{3}.Stats{'GenTq','Std'}]);
            DataCnt.GenTqSTD.std(iV,iTI,iP) = std([DataIN{1}.Stats{'GenTq','Std'},DataIN{2}.Stats{'GenTq','Std'},DataIN{3}.Stats{'GenTq','Std'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.GenTqSTD.mean(iV,iTI,iP) = 0;
                DataCnt.GenTqSTD.std(iV,iTI,iP) = 0;
            end


            DataCnt.BlPitchSTD.mean(iV,iTI,iP)  = mean([DataIN{1}.Stats{'BldPitch1','Std'},DataIN{2}.Stats{'BldPitch1','Std'},DataIN{3}.Stats{'BldPitch1','Std'}]);
            DataCnt.BlPitchSTD.std(iV,iTI,iP)  = std([DataIN{1}.Stats{'BldPitch1','Std'},DataIN{2}.Stats{'BldPitch1','Std'},DataIN{3}.Stats{'BldPitch1','Std'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.BlPitchSTD.mean(iV,iTI,iP) = 0;
                DataCnt.BlPitchSTD.std(iV,iTI,iP) = 0;
            end


            DataCnt.BlPitchTrav.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'BldPitch1','SyAbs'},DataIN{2}.Stats{'BldPitch1','SyAbs'},DataIN{3}.Stats{'BldPitch1','SyAbs'}]);
            DataCnt.BlPitchTrav.std(iV,iTI,iP) = std([DataIN{1}.Stats{'BldPitch1','SyAbs'},DataIN{2}.Stats{'BldPitch1','SyAbs'},DataIN{3}.Stats{'BldPitch1','SyAbs'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.BlPitchTrav.mean(iV,iTI,iP) = 0;
                DataCnt.BlPitchTrav.std(iV,iTI,iP) = 0;
            end


            DataCnt.TTDspFA.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'TTDspFA','Std'},DataIN{2}.Stats{'TTDspFA','Std'},DataIN{3}.Stats{'TTDspFA','Std'}]);
            DataCnt.TTDspFA.std(iV,iTI,iP) = std([DataIN{1}.Stats{'TTDspFA','Std'},DataIN{2}.Stats{'TTDspFA','Std'},DataIN{3}.Stats{'TTDspFA','Std'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.TTDspFA.mean(iV,iTI,iP) = 0;
                DataCnt.TTDspFA.std(iV,iTI,iP) = 0;
            end


            DataCnt.TBMx.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'TwrBsMxt','DEL1Hz_4'},DataIN{2}.Stats{'TwrBsMxt','DEL1Hz_4'},DataIN{3}.Stats{'TwrBsMxt','DEL1Hz_4'}]);
            DataCnt.TBMx.std(iV,iTI,iP) = std([DataIN{1}.Stats{'TwrBsMxt','DEL1Hz_4'},DataIN{2}.Stats{'TwrBsMxt','DEL1Hz_4'},DataIN{3}.Stats{'TwrBsMxt','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.TBMx.mean(iV,iTI,iP) = 0;
                DataCnt.TBMx.std(iV,iTI,iP) = 0;
            end


            DataCnt.TBMy.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'TwrBsMyt','DEL1Hz_4'},DataIN{2}.Stats{'TwrBsMyt','DEL1Hz_4'},DataIN{3}.Stats{'TwrBsMyt','DEL1Hz_4'}]);
            DataCnt.TBMy.std(iV,iTI,iP) = std([DataIN{1}.Stats{'TwrBsMyt','DEL1Hz_4'},DataIN{2}.Stats{'TwrBsMyt','DEL1Hz_4'},DataIN{3}.Stats{'TwrBsMyt','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.TBMy.mean(iV,iTI,iP) = 0;
                DataCnt.TBMy.std(iV,iTI,iP) = 0;
            end


            DataCnt.TBMz.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'TwrBsMzt','DEL1Hz_4'},DataIN{2}.Stats{'TwrBsMzt','DEL1Hz_4'},DataIN{3}.Stats{'TwrBsMzt','DEL1Hz_4'}]);
            DataCnt.TBMz.std(iV,iTI,iP) = std([DataIN{1}.Stats{'TwrBsMzt','DEL1Hz_4'},DataIN{2}.Stats{'TwrBsMzt','DEL1Hz_4'},DataIN{3}.Stats{'TwrBsMzt','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<20
                DataCnt.TBMz.mean(iV,iTI,iP) = 0;
                DataCnt.TBMz.std(iV,iTI,iP) = 0;
            end


            DataCnt.BRMx.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'RootMxb1','DEL1Hz_10'},DataIN{2}.Stats{'RootMxb1','DEL1Hz_10'},DataIN{3}.Stats{'RootMxb1','DEL1Hz_10'}]);
            DataCnt.BRMx.std(iV,iTI,iP) = std([DataIN{1}.Stats{'RootMxb1','DEL1Hz_10'},DataIN{2}.Stats{'RootMxb1','DEL1Hz_10'},DataIN{3}.Stats{'RootMxb1','DEL1Hz_10'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.BRMx.mean(iV,iTI,iP) = 0;
                DataCnt.BRMx.std(iV,iTI,iP) = 0;
            end


            DataCnt.BRMy.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'RootMyb1','DEL1Hz_10'},DataIN{2}.Stats{'RootMyb1','DEL1Hz_10'},DataIN{3}.Stats{'RootMyb1','DEL1Hz_10'}]);
            DataCnt.BRMy.std(iV,iTI,iP) = std([DataIN{1}.Stats{'RootMyb1','DEL1Hz_10'},DataIN{2}.Stats{'RootMyb1','DEL1Hz_10'},DataIN{3}.Stats{'RootMyb1','DEL1Hz_10'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.BRMy.mean(iV,iTI,iP) = 0;
                DataCnt.BRMy.std(iV,iTI,iP) = 0;
            end


            DataCnt.BRMz.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'RootMzb1','DEL1Hz_10'},DataIN{2}.Stats{'RootMzb1','DEL1Hz_10'},DataIN{3}.Stats{'RootMzb1','DEL1Hz_10'}]);
            DataCnt.BRMz.std(iV,iTI,iP) = std([DataIN{1}.Stats{'RootMzb1','DEL1Hz_10'},DataIN{2}.Stats{'RootMzb1','DEL1Hz_10'},DataIN{3}.Stats{'RootMzb1','DEL1Hz_10'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.BRMz.mean(iV,iTI,iP) = 0;
                DataCnt.BRMz.std(iV,iTI,iP) = 0;
            end


            DataCnt.BROop.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'RootMyc1','DEL1Hz_10'},DataIN{2}.Stats{'RootMyc1','DEL1Hz_10'},DataIN{3}.Stats{'RootMyc1','DEL1Hz_10'}]);
            DataCnt.BROop.std(iV,iTI,iP) = std([DataIN{1}.Stats{'RootMyc1','DEL1Hz_10'},DataIN{2}.Stats{'RootMyc1','DEL1Hz_10'},DataIN{3}.Stats{'RootMyc1','DEL1Hz_10'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.BROop.mean(iV,iTI,iP) = 0;
                DataCnt.BROop.std(iV,iTI,iP) = 0;
            end


            DataCnt.BRIp.mean(iV,iTI,iP)  = mean([DataIN{1}.Stats{'RootMxc1','DEL1Hz_10'},DataIN{2}.Stats{'RootMxc1','DEL1Hz_10'},DataIN{3}.Stats{'RootMxc1','DEL1Hz_10'}]);
            DataCnt.BRIp.std(iV,iTI,iP)  = std([DataIN{1}.Stats{'RootMxc1','DEL1Hz_10'},DataIN{2}.Stats{'RootMxc1','DEL1Hz_10'},DataIN{3}.Stats{'RootMxc1','DEL1Hz_10'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.BRIp.mean(iV,iTI,iP) = 0;
                DataCnt.BRIp.std(iV,iTI,iP) = 0;
            end


            DataCnt.TTMx.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'YawBrMxn','DEL1Hz_4'},DataIN{2}.Stats{'YawBrMxn','DEL1Hz_4'},DataIN{3}.Stats{'YawBrMxn','DEL1Hz_4'}]);
            DataCnt.TTMx.std(iV,iTI,iP) = std([DataIN{1}.Stats{'YawBrMxn','DEL1Hz_4'},DataIN{2}.Stats{'YawBrMxn','DEL1Hz_4'},DataIN{3}.Stats{'YawBrMxn','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.TTMx.mean(iV,iTI,iP) = 0;
                DataCnt.TTMx.std(iV,iTI,iP) = 0;
            end


            DataCnt.TTMy.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'YawBrMyn','DEL1Hz_4'},DataIN{2}.Stats{'YawBrMyn','DEL1Hz_4'},DataIN{3}.Stats{'YawBrMyn','DEL1Hz_4'}]);
            DataCnt.TTMy.std(iV,iTI,iP) = std([DataIN{1}.Stats{'YawBrMyn','DEL1Hz_4'},DataIN{2}.Stats{'YawBrMyn','DEL1Hz_4'},DataIN{3}.Stats{'YawBrMyn','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.TTMy.mean(iV,iTI,iP) = 0;
                DataCnt.TTMy.std(iV,iTI,iP) = 0;
            end


            DataCnt.TTMz.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'YawBrMzn','DEL1Hz_4'},DataIN{2}.Stats{'YawBrMzn','DEL1Hz_4'},DataIN{3}.Stats{'YawBrMzn','DEL1Hz_4'}]);
            DataCnt.TTMz.std(iV,iTI,iP) = std([DataIN{1}.Stats{'YawBrMzn','DEL1Hz_4'},DataIN{2}.Stats{'YawBrMzn','DEL1Hz_4'},DataIN{3}.Stats{'YawBrMzn','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.TTMz.mean(iV,iTI,iP) = 0;
                DataCnt.TTMz.std(iV,iTI,iP) = 0;
            end


            DataCnt.LSSMy.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'LSSTipMys','DEL1Hz_4'},DataIN{2}.Stats{'LSSTipMys','DEL1Hz_4'},DataIN{3}.Stats{'LSSTipMys','DEL1Hz_4'}]);
            DataCnt.LSSMy.std(iV,iTI,iP) = std([DataIN{1}.Stats{'LSSTipMys','DEL1Hz_4'},DataIN{2}.Stats{'LSSTipMys','DEL1Hz_4'},DataIN{3}.Stats{'LSSTipMys','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.LSSMy.mean(iV,iTI,iP) = 0;
                DataCnt.LSSMy.std(iV,iTI,iP) = 0;
            end


            DataCnt.LSSMz.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'LSSTipMzs','DEL1Hz_4'},DataIN{2}.Stats{'LSSTipMzs','DEL1Hz_4'},DataIN{3}.Stats{'LSSTipMzs','DEL1Hz_4'}]);
            DataCnt.LSSMz.std(iV,iTI,iP) = std([DataIN{1}.Stats{'LSSTipMzs','DEL1Hz_4'},DataIN{2}.Stats{'LSSTipMzs','DEL1Hz_4'},DataIN{3}.Stats{'LSSTipMzs','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.LSSMz.mean(iV,iTI,iP) = 0;
                DataCnt.LSSMz.std(iV,iTI,iP) = 0;
            end


            DataCnt.LSSTq.mean(iV,iTI,iP) = mean([DataIN{1}.Stats{'LSShftTq','DEL1Hz_4'},DataIN{2}.Stats{'LSShftTq','DEL1Hz_4'},DataIN{3}.Stats{'LSShftTq','DEL1Hz_4'}]);
            DataCnt.LSSTq.std(iV,iTI,iP) = std([DataIN{1}.Stats{'LSShftTq','DEL1Hz_4'},DataIN{2}.Stats{'LSShftTq','DEL1Hz_4'},DataIN{3}.Stats{'LSShftTq','DEL1Hz_4'}]);
            if  DataCnt.Power.mean(iV,iTI,iP)<Poff
                DataCnt.LSSTq.mean(iV,iTI,iP) = 0;
                DataCnt.LSSTq.std(iV,iTI,iP) = 0;
            end


            clear DataIN curfileN1 curfileN2 curfileN3 curSubFold1 curSubFold2 curSubFold3
            %             toc
        end

    end
end
DataCnt.Dimensions.Dim1{1} = 'WSP';
DataCnt.Dimensions.Dim1{2} = Vvec;
DataCnt.Dimensions.Dim2{1}  = 'TI';
DataCnt.Dimensions.Dim2{2} = TIvec;
DataCnt.Dimensions.Dim3{1}  = 'Prat';
DataCnt.Dimensions.Dim3{2} = Pvec;

save(savename,'DataCnt')
% figure,plot(Vvec,squeeze(DataCnt.Power.mean(:,2,11))),grid on
% figure,plot(Pvec,squeeze(DataCnt.Power.mean(12,2,:))),grid on
figure,plot(Pvec,squeeze(DataCnt.BRMy.mean(8,2,:))),grid on



