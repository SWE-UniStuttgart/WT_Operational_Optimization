
clearvars
close all
pathbase= 'D:\data\34_pettas\PhD\Optimization\Final\Baseline Results all\Spline\';
pathopt = 'D:\data\34_pettas\PhD\Optimization\Final\P_IBC_Shut_opt\';
OptFiles = {...
%     'DK_2017_ShutDown1' ...
%     'DK_2017_ShutDown1' ...
%     'DK_2017_IBC1' ...
%     'DK_2017_IBC2' ...
%     'DK_2017_IBC_Shut1'...
%     'DK_2017_IBC_Shut2' ...
%     'DK_2017_IBC_ShutBoosr1'...
%     'DK_all_Const80_IBC_GAmulti' ...
%     'DE_all_Const80_IBC'...
%     'DE_all_Rev_hor_all'...
% 'DK_2018_spl_ShutThres_100_IBCthres_13' ...
% 'DK_2018_spl_ShutThres_100_IBCthres_15' ...
% 'DK_2018_spl_ShutThres_100_IBCthres_17' ...
% 'DK_2018_spl_ShutThres_100_IBCthres_19' ...
% 'DK_2018_spl_ShutThres_100_IBCthres_21' ...
%          'DK_2018_spl_ShutThres_60_IBCthres_17'; ...
        %  'DK_2018_spl_ShutThres_80_IBCthres_17' ...
        %  'DK_2018_spl_ShutThres_100_IBCthres_17' ...
        %  'DK_2018_spl_ShutThres_120_IBCthres_17' ...
        %  'DK_2018_spl_ShutThres_140_IBCthres_17' ...
        %  'DK_2018_spl_ShutThres_160_IBCthres_17' ...
        %  'DK_2018_spl_ShutThres_180_IBCthres_17' ...
        %  'DK_2018_spl_ShutThres_200_IBCthres_17' ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_200' ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_250' ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_300' ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_350' ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_400' ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_450' ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_500'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_600'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_650'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_675'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_17_BoostThres_700'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_15_BoostThres_475'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_15_BoostThres_575'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_15_BoostThres_600'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_15_BoostThres_625'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_15_BoostThres_650'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_15_BoostThres_675'; ...
% 'DK_2018_spl_ShutThres_60_IBCthres_15_BoostThres_700' ...
% 'DE_2016_spl_ShutThres_20';...
% 'DE_2016_spl_ShutThres_40';...
% 'DE_2016_spl_ShutThres_60';...
% 'DE_2016_spl_ShutThres_80';...
% 'DE_2016_spl_ShutThres_100';...
% 'DE_2016_spl_ShutThres_120';...
% 'DE_2016_spl_ShutThres_140';...
% 'DE_2016_spl_ShutThres_160';...
% 'DE_2016_spl_ShutThres_180';...
% 'DE_2016_spl_ShutThres_200';...
% 'DE_2016_spl_ShutThres_220';... 
% 'DE_2016_spl_ShutThres_240';...
% 'DE_2017_spl_ShutThres_20';...
% 'DE_2017_spl_ShutThres_40';...
% 'DE_2017_spl_ShutThres_60';...
% 'DE_2017_spl_ShutThres_80';...
% 'DE_2017_spl_ShutThres_100';...
% 'DE_2017_spl_ShutThres_120';...
% 'DE_2017_spl_ShutThres_140';...
% 'DE_2017_spl_ShutThres_160';...
% 'DE_2017_spl_ShutThres_180';...
% 'DE_2017_spl_ShutThres_200';...
% 'DE_2018_spl_ShutThres_20';...
% 'DE_2018_spl_ShutThres_40';...
% 'DE_2018_spl_ShutThres_60';...
% 'DE_2018_spl_ShutThres_80';...
% 'DE_2018_spl_ShutThres_100';...
% 'DE_2018_spl_ShutThres_120';...
% 'DE_2018_spl_ShutThres_140';...
% 'DE_2018_spl_ShutThres_160';...
% 'DE_2018_spl_ShutThres_180';...
% % 'DE_2018_spl_ShutThres_200';...
% 'DE_all_spl_ShutThres_20';...
% 'DE_all_spl_ShutThres_40';...
% 'DE_all_spl_ShutThres_60';...
% 'DE_all_spl_ShutThres_80';...
% 'DE_all_spl_ShutThres_100';...
% 'DE_all_spl_ShutThres_120';...
% 'DE_all_spl_ShutThres_140';...
% 'DE_all_spl_ShutThres_160';...
% 'DE_all_spl_ShutThres_180';...
% 'DE_all_spl_ShutThres_200';...
% 'DK_2013_spl_ShutThres_20';...
% 'DK_2013_spl_ShutThres_40';...
% 'DK_2013_spl_ShutThres_60';...
% 'DK_2013_spl_ShutThres_80';...
% 'DK_2013_spl_ShutThres_100';...
% 'DK_2013_spl_ShutThres_120';...
% 'DK_2013_spl_ShutThres_140';...
% 'DK_2013_spl_ShutThres_160';...
% 'DK_2013_spl_ShutThres_180';...
% 'DK_2013_spl_ShutThres_200';...
% 'DK_2013_spl_ShutThres_220';...
% 'DK_2013_spl_ShutThres_240';...
% 'DK_2014_spl_ShutThres_20';...
% 'DK_2014_spl_ShutThres_40';...
% 'DK_2014_spl_ShutThres_60';...
% 'DK_2014_spl_ShutThres_80';...
% 'DK_2014_spl_ShutThres_100';...
% 'DK_2014_spl_ShutThres_120';...
% 'DK_2014_spl_ShutThres_140';...
% 'DK_2014_spl_ShutThres_160';...
% 'DK_2014_spl_ShutThres_180';...
% 'DK_2014_spl_ShutThres_200';...
% 'DK_2014_spl_ShutThres_220';...
% 'DK_2014_spl_ShutThres_240';...
% 'DK_2015_spl_ShutThres_20';...
% 'DK_2015_spl_ShutThres_40';...
% 'DK_2015_spl_ShutThres_60';...
% 'DK_2015_spl_ShutThres_80';...
% 'DK_2015_spl_ShutThres_100';...
% 'DK_2015_spl_ShutThres_120';...
% 'DK_2015_spl_ShutThres_140';...
% 'DK_2015_spl_ShutThres_160';...
% 'DK_2015_spl_ShutThres_180';...
% 'DK_2015_spl_ShutThres_200';...
% 'DK_2015_spl_ShutThres_220';...
% 'DK_2015_spl_ShutThres_240';...
% 'DK_2016_spl_ShutThres_20';...
% 'DK_2016_spl_ShutThres_40';...
% 'DK_2016_spl_ShutThres_60';...
% 'DK_2016_spl_ShutThres_80';...
% 'DK_2016_spl_ShutThres_100';...
% 'DK_2016_spl_ShutThres_120';...
% 'DK_2016_spl_ShutThres_140';...
% 'DK_2016_spl_ShutThres_160';...
% 'DK_2016_spl_ShutThres_180';...
% 'DK_2016_spl_ShutThres_200';...
% 'DK_2016_spl_ShutThres_220';...
% 'DK_2016_spl_ShutThres_240';...
% 'DK_2017_spl_ShutThres_20';...
% 'DK_2017_spl_ShutThres_40';...
% 'DK_2017_spl_ShutThres_60';...
% 'DK_2017_spl_ShutThres_80';...
% 'DK_2017_spl_ShutThres_100';...
% 'DK_2017_spl_ShutThres_120';...
% 'DK_2017_spl_ShutThres_140';...
% 'DK_2017_spl_ShutThres_160';...
% 'DK_2017_spl_ShutThres_180';...
% 'DK_2017_spl_ShutThres_200';...
% 'DK_2017_spl_ShutThres_220';...
% 'DK_2017_spl_ShutThres_240';...
% 'DK_2018_spl_ShutThres_20';...
% 'DK_2018_spl_ShutThres_40';...
% 'DK_2018_spl_ShutThres_60';...
% 'DK_2018_spl_ShutThres_80';...
% 'DK_2018_spl_ShutThres_100';...
% 'DK_2018_spl_ShutThres_120';...
% 'DK_2018_spl_ShutThres_140';...
% 'DK_2018_spl_ShutThres_160';...
% 'DK_2018_spl_ShutThres_180';...
% 'DK_2018_spl_ShutThres_200';...
% 'DK_2018_spl_ShutThres_220';...
% 'DK_2018_spl_ShutThres_240';...
% 'DK_2019_spl_ShutThres_20';...
% 'DK_2019_spl_ShutThres_40';...
% 'DK_2019_spl_ShutThres_60';...
% 'DK_2019_spl_ShutThres_80';...
% 'DK_2019_spl_ShutThres_100';...
% 'DK_2019_spl_ShutThres_120';...
% 'DK_2019_spl_ShutThres_140';...
% 'DK_2019_spl_ShutThres_160';...
% 'DK_2019_spl_ShutThres_180';...
% 'DK_2019_spl_ShutThres_200';...
% 'DK_2019_spl_ShutThres_220';...
% 'DK_2019_spl_ShutThres_240';...
% 'DK_2020_spl_ShutThres_20';...
% 'DK_2020_spl_ShutThres_40';...
% 'DK_2020_spl_ShutThres_60';...
% 'DK_2020_spl_ShutThres_80';...
% 'DK_2020_spl_ShutThres_100';...
% 'DK_2020_spl_ShutThres_120';...
% 'DK_2020_spl_ShutThres_140';...
% 'DK_2020_spl_ShutThres_160';...
% 'DK_2020_spl_ShutThres_180';...
% 'DK_2020_spl_ShutThres_200';...
% 'DK_2020_spl_ShutThres_220';...
% 'DK_2020_spl_ShutThres_240';...
% 'DK_all_spl_ShutThres_20';...
% 'DK_all_spl_ShutThres_40';...
% 'DK_all_spl_ShutThres_60';...
% 'DK_all_spl_ShutThres_80';...
% 'DK_all_spl_ShutThres_100';...
% 'DK_all_spl_ShutThres_120';...
% 'DK_all_spl_ShutThres_140';...
% 'DK_all_spl_ShutThres_160';...
% 'DK_all_spl_ShutThres_180';...
% 'DK_all_spl_ShutThres_200';...
% 'DK_all_spl_ShutThres_220';...
% 'DK_all_spl_ShutThres_240';...
% 'DE_2016_spl_IBCThres_11';...
% 'DE_2016_spl_IBCThres_12';...
% 'DE_2016_spl_IBCThres_13';...
% 'DE_2016_spl_IBCThres_14';...
% 'DE_2016_spl_IBCThres_15';...
% 'DE_2016_spl_IBCThres_16';...
% 'DE_2016_spl_IBCThres_17';...
% 'DE_2016_spl_IBCThres_18';...
% 'DE_2016_spl_IBCThres_19';...
% 'DE_2016_spl_IBCThres_20';...
% 'DE_2016_spl_IBCThres_21';...
% 'DE_2017_spl_IBCThres_11';...
% 'DE_2017_spl_IBCThres_12';...
% 'DE_2017_spl_IBCThres_13';...
% 'DE_2017_spl_IBCThres_14';...
% 'DE_2017_spl_IBCThres_15';...
% 'DE_2017_spl_IBCThres_16';...
% 'DE_2017_spl_IBCThres_17';...
% 'DE_2017_spl_IBCThres_18';...
% 'DE_2017_spl_IBCThres_19';...
% 'DE_2017_spl_IBCThres_20';...
% 'DE_2017_spl_IBCThres_21';...
% 'DE_2018_spl_IBCThres_11';...
% 'DE_2018_spl_IBCThres_12';...
% 'DE_2018_spl_IBCThres_13';...
% 'DE_2018_spl_IBCThres_14';...
% 'DE_2018_spl_IBCThres_15';...
% 'DE_2018_spl_IBCThres_16';...
% 'DE_2018_spl_IBCThres_17';...
% 'DE_2018_spl_IBCThres_18';...
% 'DE_2018_spl_IBCThres_19';...
% 'DE_2018_spl_IBCThres_20';...
% 'DE_2018_spl_IBCThres_21';...
% 'DE_all_spl_IBCThres_11';...
% 'DE_all_spl_IBCThres_12';...
% 'DE_all_spl_IBCThres_13';...
% 'DE_all_spl_IBCThres_14';...
% 'DE_all_spl_IBCThres_15';...
% 'DE_all_spl_IBCThres_16';...
% 'DE_all_spl_IBCThres_17';...
% 'DE_all_spl_IBCThres_18';...
% 'DE_all_spl_IBCThres_19';...
% 'DE_all_spl_IBCThres_20';...
% 'DE_all_spl_IBCThres_21';...
% 'DK_all_spl_IBCThres_11';...
% 'DK_all_spl_IBCThres_12';...
% 'DK_all_spl_IBCThres_13';...
% 'DK_all_spl_IBCThres_14';...
% 'DK_all_spl_IBCThres_15';...
% 'DK_all_spl_IBCThres_16';...
% 'DK_all_spl_IBCThres_17';...
% 'DK_all_spl_IBCThres_18';...
% 'DK_all_spl_IBCThres_19';...
% 'DK_all_spl_IBCThres_20';...
% 'DK_all_spl_IBCThres_21';...
% 'DK_2013_spl_IBCThres_11';...
% 'DK_2013_spl_IBCThres_12';...
% 'DK_2013_spl_IBCThres_13';...
% 'DK_2013_spl_IBCThres_14';...
% 'DK_2013_spl_IBCThres_15';...
% 'DK_2013_spl_IBCThres_16';...
% 'DK_2013_spl_IBCThres_17';...
% 'DK_2013_spl_IBCThres_18';...
% 'DK_2013_spl_IBCThres_19';...
% 'DK_2013_spl_IBCThres_20';...
% 'DK_2013_spl_IBCThres_21';...
% 'DK_2014_spl_IBCThres_11';...
% 'DK_2014_spl_IBCThres_12';...
% 'DK_2014_spl_IBCThres_13';...
% 'DK_2014_spl_IBCThres_14';...
% 'DK_2014_spl_IBCThres_15';...
% 'DK_2014_spl_IBCThres_16';...
% 'DK_2014_spl_IBCThres_17';...
% 'DK_2014_spl_IBCThres_18';...
% 'DK_2014_spl_IBCThres_19';...
% 'DK_2014_spl_IBCThres_20';...
% 'DK_2014_spl_IBCThres_21';...
% 'DK_2015_spl_IBCThres_11';...
% 'DK_2015_spl_IBCThres_12';...
% 'DK_2015_spl_IBCThres_13';...
% 'DK_2015_spl_IBCThres_14';...
% 'DK_2015_spl_IBCThres_15';...
% 'DK_2015_spl_IBCThres_16';...
% 'DK_2015_spl_IBCThres_17';...
% 'DK_2015_spl_IBCThres_18';...
% 'DK_2015_spl_IBCThres_19';...
% 'DK_2015_spl_IBCThres_20';...
% 'DK_2015_spl_IBCThres_21';...
% 'DK_2016_spl_IBCThres_11';...
% 'DK_2016_spl_IBCThres_12';...
% 'DK_2016_spl_IBCThres_13';...
% 'DK_2016_spl_IBCThres_14';...
% 'DK_2016_spl_IBCThres_15';...
% 'DK_2016_spl_IBCThres_16';...
% 'DK_2016_spl_IBCThres_17';...
% 'DK_2016_spl_IBCThres_18';...
% 'DK_2016_spl_IBCThres_19';...
% 'DK_2016_spl_IBCThres_20';...
% 'DK_2016_spl_IBCThres_21';...
% 'DK_2017_spl_IBCThres_11';...
% 'DK_2017_spl_IBCThres_12';...
% 'DK_2017_spl_IBCThres_13';...
% 'DK_2017_spl_IBCThres_14';...
% 'DK_2017_spl_IBCThres_15';...
% 'DK_2017_spl_IBCThres_16';...
% 'DK_2017_spl_IBCThres_17';...
% 'DK_2017_spl_IBCThres_18';...
% 'DK_2017_spl_IBCThres_19';...
% 'DK_2017_spl_IBCThres_20';...
% 'DK_2017_spl_IBCThres_21';...
% 'DK_2018_spl_IBCThres_11';...
% 'DK_2018_spl_IBCThres_12';...
% 'DK_2018_spl_IBCThres_13';...
% 'DK_2018_spl_IBCThres_14';...
% 'DK_2018_spl_IBCThres_15';...
% 'DK_2018_spl_IBCThres_16';...
% 'DK_2018_spl_IBCThres_17';...
% 'DK_2018_spl_IBCThres_18';...
% 'DK_2018_spl_IBCThres_19';...
% 'DK_2018_spl_IBCThres_20';...
% 'DK_2018_spl_IBCThres_21';...
% 'DK_2019_spl_IBCThres_11';...
% 'DK_2019_spl_IBCThres_12';...
% 'DK_2019_spl_IBCThres_13';...
% 'DK_2019_spl_IBCThres_14';...
% 'DK_2019_spl_IBCThres_15';...
% 'DK_2019_spl_IBCThres_16';...
% 'DK_2019_spl_IBCThres_17';...
% 'DK_2019_spl_IBCThres_18';...
% 'DK_2019_spl_IBCThres_19';...
% 'DK_2019_spl_IBCThres_20';...
% 'DK_2019_spl_IBCThres_21';...
% 'DK_2020_spl_IBCThres_11';...
% 'DK_2020_spl_IBCThres_12';...
% 'DK_2020_spl_IBCThres_13';...
% 'DK_2020_spl_IBCThres_14';...
% 'DK_2020_spl_IBCThres_15';...
% 'DK_2020_spl_IBCThres_16';...
% 'DK_2020_spl_IBCThres_17';...
% 'DK_2020_spl_IBCThres_18';...
% 'DK_2020_spl_IBCThres_19';...
% 'DK_2020_spl_IBCThres_20';...
% 'DK_2020_spl_IBCThres_21';...
% 'DE_2016_spl_BoostThres_200';...
% 'DE_2016_spl_BoostThres_250';...
% 'DE_2016_spl_BoostThres_300';...
% 'DE_2016_spl_BoostThres_350';...
% 'DE_2016_spl_BoostThres_400';...
% 'DE_2016_spl_BoostThres_450';...
% 'DE_2016_spl_BoostThres_500';...
% 'DE_2016_spl_BoostThres_550';...
% 'DE_2016_spl_BoostThres_600';...
% 'DE_2017_spl_BoostThres_200';...
% 'DE_2017_spl_BoostThres_250';...
% 'DE_2017_spl_BoostThres_300';...
% 'DE_2017_spl_BoostThres_350';...
% 'DE_2017_spl_BoostThres_400';...
% 'DE_2017_spl_BoostThres_450';...
% 'DE_2017_spl_BoostThres_500';...
% 'DE_2017_spl_BoostThres_550';...
% 'DE_2017_spl_BoostThres_600';...
% 'DE_2018_spl_BoostThres_200';...
% 'DE_2018_spl_BoostThres_250';...
% 'DE_2018_spl_BoostThres_300';...
% 'DE_2018_spl_BoostThres_350';...
% 'DE_2018_spl_BoostThres_400';...
% 'DE_2018_spl_BoostThres_450';...
% 'DE_2018_spl_BoostThres_500';...
% 'DE_2018_spl_BoostThres_550';...
% 'DE_2018_spl_BoostThres_600';...
% 'DE_all_spl_BoostThres_200';...
% 'DE_all_spl_BoostThres_250';...
% 'DE_all_spl_BoostThres_300';...
% 'DE_all_spl_BoostThres_350';...
% 'DE_all_spl_BoostThres_400';...
% 'DE_all_spl_BoostThres_450';...
% 'DE_all_spl_BoostThres_500';...
% 'DE_all_spl_BoostThres_550';...
% 'DE_all_spl_BoostThres_600';...
% 'DK_all_spl_BoostThres_200';...
% 'DK_all_spl_BoostThres_250';...
% 'DK_all_spl_BoostThres_300';...
% 'DK_all_spl_BoostThres_350';...
% 'DK_all_spl_BoostThres_400';...
% 'DK_all_spl_BoostThres_450';...
% 'DK_all_spl_BoostThres_500';...
% 'DK_all_spl_BoostThres_550';...
% 'DK_all_spl_BoostThres_600';...
% 'DK_2013_spl_BoostThres_200';...
% 'DK_2013_spl_BoostThres_250';...
% 'DK_2013_spl_BoostThres_300';...
% 'DK_2013_spl_BoostThres_350';...
% 'DK_2013_spl_BoostThres_400';...
% 'DK_2013_spl_BoostThres_450';...
% 'DK_2013_spl_BoostThres_500';...
% 'DK_2013_spl_BoostThres_550';...
% 'DK_2013_spl_BoostThres_600';...
% 'DK_2014_spl_BoostThres_200';...
% 'DK_2014_spl_BoostThres_250';...
% 'DK_2014_spl_BoostThres_300';...
% 'DK_2014_spl_BoostThres_350';...
% 'DK_2014_spl_BoostThres_400';...
% 'DK_2014_spl_BoostThres_450';...
% 'DK_2014_spl_BoostThres_500';...
% 'DK_2014_spl_BoostThres_550';...
% 'DK_2014_spl_BoostThres_600';...
% 'DK_2015_spl_BoostThres_200';...
% 'DK_2015_spl_BoostThres_250';...
% 'DK_2015_spl_BoostThres_300';...
% 'DK_2015_spl_BoostThres_350';...
% 'DK_2015_spl_BoostThres_400';...
% 'DK_2015_spl_BoostThres_450';...
% 'DK_2015_spl_BoostThres_500';...
% 'DK_2015_spl_BoostThres_550';...
% 'DK_2015_spl_BoostThres_600';...
% 'DK_2016_spl_BoostThres_200';...
% 'DK_2016_spl_BoostThres_250';...
% 'DK_2016_spl_BoostThres_300';...
% 'DK_2016_spl_BoostThres_350';...
% 'DK_2016_spl_BoostThres_400';...
% 'DK_2016_spl_BoostThres_450';...
% 'DK_2016_spl_BoostThres_500';...
% 'DK_2016_spl_BoostThres_550';...
% 'DK_2016_spl_BoostThres_600';...
% 'DK_2017_spl_BoostThres_200';...
% 'DK_2017_spl_BoostThres_250';...
% 'DK_2017_spl_BoostThres_300';...
% 'DK_2017_spl_BoostThres_350';...
% 'DK_2017_spl_BoostThres_400';...
% 'DK_2017_spl_BoostThres_450';...
% 'DK_2017_spl_BoostThres_500';...
% 'DK_2017_spl_BoostThres_550';...
% 'DK_2017_spl_BoostThres_600';...
% 'DK_2018_spl_BoostThres_200';...
% 'DK_2018_spl_BoostThres_250';...
% 'DK_2018_spl_BoostThres_300';...
% 'DK_2018_spl_BoostThres_350';...
% 'DK_2018_spl_BoostThres_400';...
% 'DK_2018_spl_BoostThres_450';...
% 'DK_2018_spl_BoostThres_500';...
% 'DK_2018_spl_BoostThres_550';...
% 'DK_2018_spl_BoostThres_600';...
% 'DK_2019_spl_BoostThres_200';...
% 'DK_2019_spl_BoostThres_250';...
% 'DK_2019_spl_BoostThres_300';...
% 'DK_2019_spl_BoostThres_350';...
% 'DK_2019_spl_BoostThres_400';...
% 'DK_2019_spl_BoostThres_450';...
% 'DK_2019_spl_BoostThres_500';...
% 'DK_2019_spl_BoostThres_550';...
% 'DK_2019_spl_BoostThres_600';...
% 'DK_2020_spl_BoostThres_200';...
% 'DK_2020_spl_BoostThres_250';...
% 'DK_2020_spl_BoostThres_300';...
% 'DK_2020_spl_BoostThres_350';...
% 'DK_2020_spl_BoostThres_400';...
% 'DK_2020_spl_BoostThres_450';...
% 'DK_2020_spl_BoostThres_500';...
% 'DK_2020_spl_BoostThres_550';...
% 'DK_2020_spl_BoostThres_600';...
% 'DE_2016_spl_BoostThres12_200';...
% 'DE_2016_spl_BoostThres12_250';...
% 'DE_2016_spl_BoostThres12_300';...
% 'DE_2016_spl_BoostThres12_350';...
% 'DE_2016_spl_BoostThres12_400';...
% 'DE_2016_spl_BoostThres12_450';...
% 'DE_2016_spl_BoostThres12_500';...
% 'DE_2016_spl_BoostThres12_550';...
% 'DE_2016_spl_BoostThres12_600';...
% 'DE_2016_spl_BoostThres11_400';...
% 'DE_2016_spl_BoostThres11_450';...
% 'DE_2016_spl_BoostThres13_350';...
% 'DE_2016_spl_BoostThres13_400';...
% 'DE_2016_spl_BoostThres13_450';...
% 
% 'DE_2016_spl_BoostThres13_150';...
% 'DE_2016_spl_BoostThres13_175';...
% 'DE_2016_spl_BoostThres13_200';...
% 'DE_2016_spl_BoostThres13_225';...
% 'DE_2016_spl_BoostThres13_250';...
% 
% 'DE_2016_spl_BoostThres13_275';...
% 'DE_2016_spl_BoostThres13_300';...
% 'DE_2016_spl_BoostThres13_325';...
% 'DE_2016_spl_BoostThres13_350';...
% 'DE_2016_spl_BoostThres13_375';...
% 'DE_2016_spl_BoostThres13_400';...
% 'DE_2016_spl_BoostThresHigh_200_BoostThresMid_250_Pmid_10.5';...
% 'DE_2016_spl_Ptrans_10.5_BoostThresHigh_200_BoostThresMid_200_Pmid_10.5';...
% 'DE_2016_spl_BoostThres_300_275_Good';...
% 'DE_2016_spl_BoostThres_290_200_14_11.5';...
% 'DE_2016_spl_BoostShutBC_Thres_190_130_13_11.5_14_Good';...
% 'DE_2017_spl_ShutIBCBoost_LoadNeutral-3.5REV';...
% 'DE_2017_spl_ShutIBCBoost_LoadNeutral-5REV';...
% 'DE_2017_spl_ShutIBCBoost_MaxLoadRed-0REV';...
% 'DE_2017_spl_ShutIBCBoost_LoadNeutral_MInIBCact-3.6REV';...
% 'DE_2017_spl_ShutIBCBoost_LoadCap5-8REV';...
'DE_2016_spl_ShutIBCBoost_LoadNeutral-4REV';...
% 'DE_2016_spl_ShutIBCBoost_MaxLoadRed-0REV';...
'DE_2016_spl_ShutIBCBoost_LoadNeutral_MInIBCact-3.4REV';...
'DE_2016_spl_ShutIBCBoost_LoadCap5-7REV';...
    }';

%%%%% Dirty patch to get all
clear OptFiles
% OptiflesAll=dir(fullfile(pathopt,'*DE_2016*mid*.mat'));
OptiflesAll=dir(fullfile(pathopt,'*DK_2016*.mat'));
for i =1:length(OptiflesAll)
    OptFiles{i,1}= OptiflesAll(i).name;
end


BaseFile = 'Baseline_DK_2016_spl';
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
if Opt{1,1}.Output.options.const_price.activate  ==1
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

