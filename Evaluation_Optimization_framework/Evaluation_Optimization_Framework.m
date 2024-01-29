%
%
% Evaluation of accumulated response of the wind turbine over a time series
% of wind conditions and prices. For each time step of the series a control
% mode is provided (through optimization or as a constant). The 
% response of the turbine for this control mode in terms of power, loads and actuator metrics is
% then added to the accumulated values until the end of the time series.
%
%
% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart
%
%
clearvars
close all
clc

fullpathHere = mfilename('fullpath');
pathHere = fullpathHere(1:end-34);
cd(pathHere)
addpath([pwd '\Surrogate_models'])
addpath([pwd '\Objective_functions'])
addpath([pwd '\Functions'])
addpath([pwd '\Surrogate_models'])
addpath([pwd '\Wind_Price_data'])
addpath([pwd '\Results'])

%% INPUTS

%%% Paths and files
% Surrogates
options.files.surrogate.path              = [pwd '\Surrogate_models']; % general path to surrogate models
options.files.surrogate.extensionPoly     = '\constTSR_smooth'; % 3D matrix for spline -plolynomial regression (includes also the IPC)
options.files.surrogate.extensionPoly_IBC = '\constTSR_IBC_smooth'; % 3D matrix for spline -plolynomial regression (includes also the IPC)
options.files.surrogate.extensionGPR      = '\'; % GPR surrogate models per sensor no IBC
options.files.surrogate.extensionGPR_IBC  = '\'; % GPR surrogate models per sensor with IBC
options.files.surrogate.case              = 1; % 1:spline / 2:GPR
% Data time series
options.files.TS_data  = [pwd '\Wind_Price_data\DE_all.mat']; % V m/s,TI %,Price eur/MWh,Time in hours as separate variables in columns
% Baseline result from this code for comparison
options.files.BaselineResult = [pwd '\Results\Baseline_DE_all_spl_const80.mat']; % baseline matfile to compare (if exists)6
% Flags
flags.CompareToBaseline = 1; %Switch on to prdict metrics in comparison to the defined baseline
flags.save.on = 1;  % if 1 save all the accumulated values
flags.save.savematname = [pwd '\Results\Test.mat'];% if savemat flag is on

%%% Settings control and optimization
%flag to choose between a constant control mode or flexible
options.control.mode = 10;
% 1:constant all no optimization choose power level and IBC activation no optimizer
% 2:Based on selective shut down only and manual thresholds no optimizer
% 3:Based on IBC only and manual thresholds no optimizer
% 4:Shut down+IBC manual thresholds no optimizer
% 5:Power level only based on thresholds with forecasts no optimizer
% 6:Power level+shutdown+IBC optimization with forecasts
% 7:Power level+shutdown+IBC with horizon and past information !!! not working yet !!!!
% 8:Power level+shutdown+IBC with manual thresholds no optimizer
% 9:Power level+shutdown+IBC optimization with the distribution based approach
% 10:Applying the result from 9 to fluctuating prices

% Fluctuating or constant prices
options.control.const_price.activate = 1; % If on, all prices will be replaced with a constant value to emulate subsidized operation
options.control.const_price.val = 80; % euro/MWh  constant value to emulate subsidized operation

% Set the windspeed above which the IBC loop is activated when IBC is generally on. Relevant only for tracking the IBC activation
options.control.IBC_generalV = 11; % m/s

%% Control modes options

% 1:constant
options.control.case{1}.Prat = 12;  % Rated power for constant mode (5:13 for now)
options.control.case{1}.IBC  = 0;   % 0 or 1 IBC on off

% 2:shut down only
options.control.case{2}.Prat  = 10;  % Rated power for constant mode (5:13 for now)
options.control.case{2}.IBC   = 0;   % 0 or 1 IPC on off
options.control.case{2}.logic.shut = 'Shut_down_opt';   % function that contains the thresholds and logic for shut down. Inputs: the current conditions, Output: shut down or not
options.control.case{2}.FunSettings.revThresLowTI_Shut = 60; % eur projected revenue to shut down in low TI
options.control.case{2}.FunSettings.revThresHighTI_Shut = 70; % eur  projected revenue to shut down in higher TI
options.control.case{2}.FunSettings.TI_Thres_Shut = 7; % TI threshold in %

%3:IBC only
options.control.case{3}.Prat  = 10;  % Prating for constant mode (5:13 for now)
options.control.case{3}.logic.IBC = 'IBC_opt';   % function that contains the thresholds and logic for activating IBC. Inputs: the current conditions, Output: IPC on or off
options.control.case{3}.FunSettings.VThresLowTI_IBC = 18; % m/s wind speed threshold to turn on IBC in low TI
options.control.case{3}.FunSettings.VThresHighTI_IBC = 16; % m/s wind speed threshold to turn on IBC in higher TI
options.control.case{3}.FunSettings.TI_Thres_IBC = 7; % TI threshold in %

% 4:shut down+IBC
options.control.case{4}.Prat  = 10;  % Prating for constant mode (5:13 for now)
options.control.case{4}.logic.shut = 'Shut_down_opt';   % function that contains the thresholds and logic for shut down. Inputs: the current conditions, Output: shut down or not
options.control.case{4}.logic.IBC  = 'IBC_opt';   % function that contains the thresholds and logic for activating IBC. Inputs: the current conditions, Output: IPC on or off
options.control.case{4}.FunSettings.revThresLowTI_Shut = 60; % eur
options.control.case{4}.FunSettings.revThresHighTI_Shut = 70; % eur
options.control.case{4}.FunSettings.TI_Thres_Shut = 7; % TI in %
options.control.case{4}.FunSettings.VThresLowTI_IBC = 18; % m/s
options.control.case{4}.FunSettings.VThresHighTI_IBC = 16; % m/s
options.control.case{4}.FunSettings.TI_Thres_IBC = 7; % TI in %

% 5:Optimize power level only with horizon
options.control.case{5}.logic.Prat = 'Optimize_P_Sorting_MultiSensor';   % function that contains the thresholds and logic for choosing power rating for a given horizon. Inputs: the current conditions, Output: power level
options.control.case{5}.horizon_duration = 'all';%7*24; % forecast horizon length in hours or all to make one block of the whole time series
options.control.case{5}.logic.bins  = 1; % if 0 the code expects optimization of each V-price bin value if 0 it expects power ratings for each time step of the block
options.control.case{5}.logic.baseP = 10; % 5-13 chooses the power rating for the baseline case
options.control.case{5}.logic.VbinSize  = 2; % m/s size for V bin discretization
options.control.case{5}.logic.PrbinSize = 5; % eur size for Price bin discretization
options.control.case{5}.IBC = 0; % switch IBC on
options.control.case{5}.FunSettings.percToTotal = 0; % The percentage of contribution of the most contributing bins per objective. If 0 no sorting
options.control.case{5}.FunSettings.wDAM = 0.2; %weigthing in objective function for DEL delta
options.control.case{5}.FunSettings.wREV = 0.8; %weigthing in objective function for revenue delta
options.control.case{5}.FunSettings.PenaltyRev = 0.5; % percentage value of reduction in revenue above which penalty is applied to the output of the objective function
options.control.case{5}.FunSettings.PenaltyDam = 5;%  percentage value of increase in revenue above which penalty is applied to the output of the objective function
options.control.case{5}.FunSettings.loadsens = {'TBMx' 'BRMy' 'TBMy' 'LSSTq' 'BRMz' 'TBMy' };%{'BROop' 'LSSTq' 'TBMy' 'TBMx' 'BRMz' 'TTMy'};
options.control.case{5}.FunSettings.Method = 'GA';  % choose optimizer 'fmincon'  'GA'  'GAmulti' 'GlobalSearch' 'MultiStart'

% 6:Rating+shutdown+IBC with horizon --The most general case for forecasts
options.control.case{6}.logic.Prat = 'Optimize_P_Shut_IBC_Sorting_MultiSensor';   % function that contains the thresholds and logic for choosing power rating for a given horizon. Inputs: the current conditions, Output: power level
options.control.case{6}.horizon_duration = 24; % 'all'   % forecast horizon length in hours <= total length of TS or 'all'
options.control.case{6}.logic.bins  = 0; % if 0 the code expects optimization of each value of the time series if 1 it expects power ratings per bin
options.control.case{6}.logic.baseP = 10; % 5-13 chooses the power rating for the baseline case
options.control.case{6}.logic.VbinSize = 1; % m/s size for V bin discretization
options.control.case{6}.logic.PrbinSize = 5; % eur size for Price bin discretization
options.control.case{6}.FunSettings.percToTotal = 0.8; % [0 1]The percentage of contribution of the most contributing bins per objective. If 0 no sorting
options.control.case{6}.FunSettings.wDAM = 0.2; %weigthing in objective function for DEL delta
options.control.case{6}.FunSettings.wREV = 0.8; %weigthing in objective function for revenue delta
options.control.case{6}.FunSettings.dlim = -5; %maximum damage reduction allowed in the objective function. Heuristic to push the optimizer away from large load reductions
options.control.case{6}.FunSettings.revlim = 30; % maximum revenue increase allowed in the objective function. Heuristic to push the optimizer away from large revenue increases
options.control.case{6}.FunSettings.PenaltyRev = 0.2; % percentage value of increase/reduction in revenue above which penalty is applied to the output of the objective function
options.control.case{6}.FunSettings.PenaltyDam = 1;%  percentage value of increase/reduction in loads above which penalty is applied to the output of the objective function
options.control.case{6}.FunSettings.loadsens = { 'BRMx' 'BRMy' 'BROop' 'BRIp' 'TTMy' 'BRMy'  'TBMy' 'BRMz' 'BRMz' 'LSSTq' 'TBMx'};%{'BROop' 'LSSTq' 'TBMy' 'TBMx' 'BRMz' 'TTMy'};
options.control.case{6}.FunSettings.Method = 'ParticleSwarm';  % choose optimizer 'fmincon'  'GA'  'GAmulti' 'GlobalSearch' 'PatternSearch' 'ParticleSwarm'
options.control.case{6}.FunSettings.MaxP = 13; % maximum level of boosting (10 13]
options.control.case{6}.FunSettings.constval.Shut = 1; % 0 optimize for shut down threshold / 1 use a constant value
options.control.case{6}.FunSettings.constval.ShutVal = 0; % revenue value to shut if the previous flag is on
options.control.case{6}.FunSettings.constval.IBC = 1; % 0 optimize for IBC threshold / 1 use a constant value
options.control.case{6}.FunSettings.constval.IBCVal = 34; % wsp value to turn IBC on if the previous flag is on

% 7: Rating+shutdown+IBC with horizon and past information -!- Under Construction / Not usable at the moment-!-
options.control.case{7}.logic.Prat = 'Power_rating_opt';  % function that contains the thresholds and logic for choosing power rating for a given horizon. Inputs: the current conditions, Output: power level
options.control.case{7}.logic.shut = 'Shut_down_opt';    % function that contains the thresholds and logic for shut down. Inputs: the current conditions, Output: shut down or not
options.control.case{7}.logic.IBC  = 'IBC_opt';   % function that contains the thresholds and logic for activating IBC. Inputs: the current conditions, Output: IPC on or off
options.control.case{7}.horizon_duration = 30*24;   % forecast horizon length in hours <= total length of TS
options.control.case{7}.past_duration    = 365*24;   % duration of past window in hours >= 0
options.control.case{7}.logic.bins  = 1;  % if 0 the code expects optimization of each value of the time series if 0 it expects power ratings for each time step of the block
options.control.case{7}.logic.baseP = 10; % 5-13 chooses the power rating for the baseline case
options.control.case{7}.logic.VbinSize  = 0.5; % m/s size for V bin discretization
options.control.case{7}.logic.PrbinSize = 1;   % eur size for Price bin discretization

% 8:shut down +IBC+ Boost [manual optimization/no optimizer]
options.control.case{8}.Prat  = 10;  % Prating for constant mode (5:13 for now)
options.control.case{8}.logic.shut = 'Shut_down_opt';   % function that contains the thresholds and logic for shut down. Inputs: the current conditions, Output: shut down or not
options.control.case{8}.logic.IBC  = 'IBC_opt';   % function that contains the thresholds and logic for activating IBC. Inputs: the current conditions, Output: IPC on or off
options.control.case{8}.logic.Boost = 'Boost_opt';   % function that contains the thresholds and logic for turning on power boosting. Inputs: the current conditions, Output: Power boost yes/no what level
options.control.case{8}.FunSettings.revThresLowTI_Shut = 30; % eur shut thres for low ti
options.control.case{8}.FunSettings.revThresHighTI_Shut = 25; % eur  eur shut thres for high ti
options.control.case{8}.FunSettings.TI_Thres_Shut = 6; % % definition of high ti for shut
options.control.case{8}.FunSettings.VThresLowTI_IBC = 19.5; % m/s IBC act threshold for low TI
options.control.case{8}.FunSettings.VThresHighTI_IBC = 17; % m/s IBC act threshold for high TI
options.control.case{8}.FunSettings.TI_Thres_IBC = 6; % % definition of high ti for IBC
options.control.case{8}.FunSettings.revThresHigh_Boost = 195; % eur revenue threshold to just go full on
options.control.case{8}.FunSettings.revThresHigh2_Boost = 275; % eur revenue threshold to boost
options.control.case{8}.FunSettings.Pmax = 11.4; % max power level allowed for the full boost in high wsp
options.control.case{8}.FunSettings.Pmax2 = 13; % max power level allowed for the full boost in high wsp
options.control.case{8}.FunSettings.revThresMid_Boost = 100; % eur lower threshold for smaller boosting. A ramp or other function could be also applied here
options.control.case{8}.FunSettings.Pmid = 11.5; % power level target for the lower threshold boost
options.control.case{8}.FunSettings.TI_Thres_Boost = 3.0; % % definition of low ti for boosting (not used for now)
options.control.case{8}.FunSettings.VThresHigh_Boost = 14.6; % wsp threshold after which boosting applies
options.control.case{8}.FunSettings.VThresHigh2_Boost = 17; % wsp second threshold after which second boosting applies
options.control.case{8}.FunSettings.Ptrans = 10; % power level for low above rated wind speeds below the previous thrshold. Here loads are very senitive better stay at 10
options.control.case{8}.FunSettings.Plow = 9.7 ; % pwr level for wind speeds below the wsp and revenue boost thresholds. Can be used to redule load in relativelt lower revenue

% 9: P, Boost, Shut, IPC based on wsp probability (e.g weibull) for CONSTANT prices
options.control.case{9}.logic.Prat = 'Optimize_P_Shut_IBC_ConstantPrices_MultiSensor_WeibOnly';
options.control.case{9}.logic.baseP = 10; % 5-13 chooses the power rating for the baseline case
options.control.case{9}.logic.VbinSize = 0.5; % m/s size for V bin discretization
options.control.case{9}.data.Weibull = [];%[11.8 2.79] % provide the long term weibul pdf/cdf (in terms of [ScaleParameter ShapeParameter]) from long term measurements (e.g mean of all DK years). If empty it will be  created by the data
options.control.case{9}.IBC = 0; % switch IBC on
options.control.case{9}.FunSettings.wDAM = 0.1; %weigthing in objective function for DEL delta
options.control.case{9}.FunSettings.wREV = 0.9; %weigthing in objective function for revenue delta
options.control.case{9}.FunSettings.dlim = -5; %maximum damage reduction allowed in the objective function. Heuristic to push the optimizer away from large load reductions
options.control.case{9}.FunSettings.revlim = 20; % maximum revenue increase allowed in the objective function. Heuristic to push the optimizer away from large revenue increases
options.control.case{9}.FunSettings.PenaltyRev = 0.1; % percentage value of reduction in revenue above which penalty is applied to the output of the objective function
options.control.case{9}.FunSettings.PenaltyDam = 1;%  percentage value of increase in damage above which penalty is applied to the output of the objective function
options.control.case{9}.FunSettings.MaxP   = 11; % maximum level of boosting [5 13]
options.control.case{9}.FunSettings.ConsTI   = 8; % constant level of ti to be considered for the calculations
options.control.case{9}.FunSettings.loadsens = { 'TBMy' 'BROop' 'BROop' 'BRMz' 'BRMy' 'BRMy' 'TTMy' 'BRIp'  'BRMz'  };%{'BROop' 'LSSTq' 'TBMy' 'TBMx' 'BRMz' 'TTMy'};
options.control.case{9}.FunSettings.Method = 'GA';  % choose optimizer 'fmincon'  'GA'  'GAmulti' 'GlobalSearch' 'ParticleSwarm'
options.control.case{9}.FunSettings.constval.Shut = 0; % 0 optimize for shut down threshold / 1 use a constant value
options.control.case{9}.FunSettings.constval.ShutVal = 0; % revenue value to shut if the previous flag is on
options.control.case{9}.FunSettings.constval.IBC = 0; % 0 optimize for IBC threshold / 1 use a constant value
options.control.case{9}.FunSettings.constval.IBCVal = 11; % wsp value to turn IBC on if the previous flag is on

% 10: Testimg the optimal points based on wsp probability (e.g weibull) for
% CONSTANT prices in fluctuating prices.
options.control.case{10}.logic.Prat = 'Optimize_P_Shut_IBC_WeibOnly';
options.control.case{10}.logic.baseP = 10; % 5-13 chooses the power rating for the baseline case
options.control.case{10}.logic.VbinSize = 0.5; % m/s size for V bin discretization
options.control.case{10}.IBC = 0; % switch IBC on
options.control.case{10}.FunSettings.ResFile = [pwd '\Results\DE_All_maxLoadRed_P130_Weib_constTSR.mat']; % result file full path

%% Load stuff and create variables

% Load surrogate models
if options.files.surrogate.case==1   % Spline
    Surrogate{1} = load([options.files.surrogate.path options.files.surrogate.extensionPoly]);     % no IPC
    Surrogate{2} = load([options.files.surrogate.path options.files.surrogate.extensionPoly_IBC]); % yes IPC
elseif options.files.surrogate.case==2   % GPR
    Surrogate{1} = load([options.files.surrogate.path options.files.surrogate.extensionGPR]); 
    Surrogate{2} = load([options.files.surrogate.path options.files.surrogate.extensionGPR_IBC]);
end

% Load data time series
Data_TS = load(options.files.TS_data);
if options.control.const_price.activate == 1
    Data_TS.Price(:) = options.control.const_price.val;
end
if options.control.mode == 11
    Data_TS_unc = load(options.control.case{11}.Unc_DataSet);
end

%% ugly clean up

%Ugly clean ups, should be done before feeding here. But hey, shit happens!
if ~isempty(find(isnan(Data_TS.Price))) && length(find(isnan(Data_TS.Price)))<4
    Data_TS.Price(isnan(Data_TS.Price)) = interp1(Data_TS.Time{1,2}(~isnan(Data_TS.Price)),Data_TS.Price(~isnan(Data_TS.Price)),Data_TS.Time{1,2}(isnan(Data_TS.Price)));
elseif  ~isempty(find(isnan(Data_TS.Price))) %#ok<*EFIND>
    Data_TS.Price(isnan(Data_TS.Price)) = 0;
end
if ~isempty(find(isnan(Data_TS.V))) && length(find(isnan(Data_TS.V)))<4
    Data_TS.V(isnan(Data_TS.V)) = interp1(Data_TS.Time{1,2}(~isnan(Data_TS.V)),Data_TS.V(~isnan(Data_TS.V)),Data_TS.Time{1,2}(isnan(Data_TS.V)));
elseif  ~isempty(find(isnan(Data_TS.V))) %#ok<*EFIND>
    Data_TS.V(isnan(Data_TS.V)) = 2;
end
if ~isempty(find(isnan(Data_TS.TI))) && length(find(isnan(Data_TS.TI)))<4
    Data_TS.TI(isnan(Data_TS.TI)) = interp1(Data_TS.Time{1,2}(~isnan(Data_TS.TI)),Data_TS.V(~isnan(Data_TS.TI)),Data_TS.Time{1,2}(isnan(Data_TS.TI)));
elseif  ~isempty(find(isnan(Data_TS.TI))) %#ok<*EFIND>
    Data_TS.TI(isnan(Data_TS.TI)) = 5;
end

if options.control.mode == 11
    if ~isempty(find(isnan(Data_TS_unc.Price))) && length(find(isnan(Data_TS_unc.Price)))<4
        Data_TS_unc.Price(isnan(Data_TS_unc.Price)) = interp1(Data_TS_unc.Time{1,2}(~isnan(Data_TS_unc.Price)),Data_TS_unc.Price(~isnan(Data_TS_unc.Price)),Data_TS_unc.Time{1,2}(isnan(Data_TS_unc.Price)));
    elseif  ~isempty(find(isnan(Data_TS_unc.Price))) %#ok<*EFIND>
        Data_TS_unc.Price(isnan(Data_TS_unc.Price)) = 0;
    end
    if ~isempty(find(isnan(Data_TS_unc.V))) && length(find(isnan(Data_TS_unc.V)))<4
        Data_TS_unc.V(isnan(Data_TS_unc.V)) = interp1(Data_TS_unc.Time{1,2}(~isnan(Data_TS_unc.V)),Data_TS_unc.V(~isnan(Data_TS_unc.V)),Data_TS_unc.Time{1,2}(isnan(Data_TS_unc.V)));
    elseif  ~isempty(find(isnan(Data_TS_unc.V))) %#ok<*EFIND>
        Data_TS_unc.V(isnan(Data_TS_unc.V)) = 2;
    end
    if ~isempty(find(isnan(Data_TS_unc.TI))) && length(find(isnan(Data_TS_unc.TI)))<4
        Data_TS_unc.TI(isnan(Data_TS_unc.TI)) = interp1(Data_TS_unc.Time{1,2}(~isnan(Data_TS_unc.TI)),Data_TS_unc.V(~isnan(Data_TS_unc.TI)),Data_TS_unc.Time{1,2}(isnan(Data_TS_unc.TI)));
    elseif  ~isempty(find(isnan(Data_TS_unc.TI))) %#ok<*EFIND>
        Data_TS_unc.TI(isnan(Data_TS_unc.TI)) = 5;
    end
end

%% Split time to blocks according to the interval requested
if options.control.mode == 1 || options.control.mode == 2 || options.control.mode == 3 || options.control.mode == 4 || options.control.mode == 8 ||  options.control.mode == 9 ||  options.control.mode == 10  || ischar(options.control.case{options.control.mode}.horizon_duration)
    Data_TS_Block{1}.Price = Data_TS.Price; %#ok<*SAGROW>
    Data_TS_Block{1}.V     = Data_TS.V;
    Data_TS_Block{1}.TI    = Data_TS.TI;
    Data_TS_Block{1}.Time{1,1} = Data_TS.Time{1,1};
    Data_TS_Block{1}.Time{1,2} = Data_TS.Time{1,2};
elseif options.control.mode == 5 || options.control.mode == 6 || options.control.mode == 7 || options.control.mode == 11
    % split the data to horizons (all is just the case of 1 horizon with all the time )
    int.x = floor(length(Data_TS.Time{1,2})/(options.control.case{options.control.mode}.horizon_duration));
    for i = 1:int.x
        int.t_block = [(i-1)*options.control.case{options.control.mode}.horizon_duration+1 i*options.control.case{options.control.mode}.horizon_duration  ];
        Data_TS_Block{i}.Price = Data_TS.Price(int.t_block(1):int.t_block(2)); %#ok<*SAGROW>
        Data_TS_Block{i}.V     = Data_TS.V(int.t_block(1):int.t_block(2));
        Data_TS_Block{i}.TI    = Data_TS.TI(int.t_block(1):int.t_block(2));
        Data_TS_Block{i}.Time{1,1} = Data_TS.Time{1,1}(int.t_block(1):int.t_block(2));
        Data_TS_Block{i}.Time{1,2} = Data_TS.Time{1,2}(int.t_block(1):int.t_block(2));
        if options.control.mode == 11
            Data_TS_Block_Unc{i}.Price = Data_TS_unc.Price(int.t_block(1):int.t_block(2)); %#ok<*SAGROW>
            Data_TS_Block_Unc{i}.V     = Data_TS_unc.V(int.t_block(1):int.t_block(2));
            Data_TS_Block_Unc{i}.TI    = Data_TS_unc.TI(int.t_block(1):int.t_block(2));
            Data_TS_Block_Unc{i}.Time{1,1} = Data_TS_unc.Time{1,1}(int.t_block(1):int.t_block(2));
            Data_TS_Block_Unc{i}.Time{1,2} = Data_TS_unc.Time{1,2}(int.t_block(1):int.t_block(2));
        end
    end
    % if there is more time left at the end add it to the last block
    if int.x*options.control.case{options.control.mode}.horizon_duration<length(Data_TS.Time{1,2})
        Data_TS_Block{int.x}.Price = [Data_TS_Block{int.x}.Price;Data_TS.Price(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
        Data_TS_Block{int.x}.V     = [Data_TS_Block{int.x}.V;Data_TS.V(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
        Data_TS_Block{int.x}.TI    = [Data_TS_Block{int.x}.TI;Data_TS.TI(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
        Data_TS_Block{int.x}.Time{1,1} = [Data_TS_Block{int.x}.Time{1,1};Data_TS.Time{1,1}(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
        Data_TS_Block{int.x}.Time{1,2} = [Data_TS_Block{int.x}.Time{1,2};Data_TS.Time{1,2}(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
        if options.control.mode == 11
            Data_TS_Block_Unc{int.x}.Price = [Data_TS_Block_Unc{int.x}.Price;Data_TS_unc.Price(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
            Data_TS_Block_Unc{int.x}.V     = [Data_TS_Block_Unc{int.x}.V;Data_TS_unc.V(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
            Data_TS_Block_Unc{int.x}.TI    = [Data_TS_Block_Unc{int.x}.TI;Data_TS_unc.TI(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
            Data_TS_Block_Unc{int.x}.Time{1,1} = [Data_TS_Block_Unc{int.x}.Time{1,1};Data_TS_unc.Time{1,1}(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
            Data_TS_Block_Unc{int.x}.Time{1,2} = [Data_TS_Block_Unc{int.x}.Time{1,2};Data_TS_unc.Time{1,2}(int.x*options.control.case{options.control.mode}.horizon_duration+1:end)];
        end
    end
end
clear int
% here uncetrainty could be added to the perfect forecasts

%% Create power curves to speed up revenue calculations
if options.control.mode == 2 || options.control.mode == 4 || options.control.mode == 8
    % create the power curve and load curve to feed to the function for quick evaluations of projected revenue and damagbe
    Vvec = [4:0.25:25]';
    if options.files.surrogate.case == 1
        if options.control.case{2}.IBC == 0
            intOut = Get_Points_From_Surrogate_noLoad(Vvec,ones(length(Vvec),1)*8,ones(length(Vvec),1).*options.control.case{2}.Prat,Surrogate{1,1},'spline') ;
        else
            intOut = Get_Points_From_Surrogate_noLoad(Vvec,ones(length(Vvec),1)*8,ones(length(Vvec),1).*options.control.case{2}.Prat,Surrogate{1,2},'spline') ;
        end
    elseif options.files.surrogate.case  == 2
        if options.control.case{2}.IBC == 0
            intOut = Get_Points_From_Surrogate_GPR(Vvec,ones(length(Vvec),1)*8,ones(length(Vvec),1).*options.control.case{2}.Prat,Surrogate{1,1});
        else
            intOut  = Get_Points_From_Surrogate_GPR(Vvec,ones(length(Vvec),1)*8,ones(length(Vvec),1).*options.control.case{2}.Prat,Surrogate{1,2});
        end
    end
    PLcurve.power = [Vvec,intOut.Power.mean/1000];
    PLcurve.load  = [Vvec,intOut.BROop.mean]; %#ok<*NBRAK>
    clear intOut Vvec
elseif options.control.mode  == 6  || options.control.mode == 9  || options.control.mode == 10  || options.control.mode == 11
    % Get PL curves for multipe power levels
    Pvec = 5:0.5:13;
    Vvec = [4:0.25:25]';
    if options.files.surrogate.case == 1
        [int.X1,int.Y1,int.Z1] = ndgrid(Surrogate{1,1}.DataCnt.Dimensions.Dim1{2},Surrogate{1, 1}.DataCnt.Dimensions.Dim2{2},Surrogate{1,1}.DataCnt.Dimensions.Dim3{2});
        intOut = squeeze(interpn(int.X1,int.Y1,int.Z1,Surrogate{1,1}.DataCnt.Power.mean,repmat(Vvec,length(Pvec),1),ones(length(Vvec)*length(Pvec),1)*8,repelem(Pvec,length(Vvec))','spline'));
        intOut = reshape(intOut,[length(Vvec),length(Pvec)]);
    elseif options.files.surrogate.case  == 2
        intOut = predict(Surrogate{1,1}.Output.Power,[[repmat(Vvec,length(Pvec),1),ones(length(Vvec)*length(Pvec),1)*8,repelem(Pvec,length(Vvec))']]);
        intOut = reshape(intOut,[length(Vvec),length(Pvec)]);
    end
    PLcurve.power = intOut; % rows: Power/cols: wsp
    PLcurve.Vvec = Vvec;
    PLcurve.Pvec = Pvec;
    clear intOut Vvec Pvec int
end
%% Loop over the time blocks

time_cnt = 0; % keeping total time
% looping blocks
for ib = 1:size(Data_TS_Block,2)
    curData = Data_TS_Block{1, ib} ;
    if options.control.mode == 11
        curData_unc = Data_TS_Block_Unc{1, ib} ;
    end

    if options.control.mode == 1 ||  options.control.mode == 3
        %    - if 1,3 just continue to the next loop to go point by point or add something more if needed here

    elseif options.control.mode == 2 || options.control.mode == 4 || options.control.mode == 8
        %    - if 2 4 8 just continue to the next loop to go point by point or add something more if needed here

    elseif options.control.mode  == 5
        %    - feed the block to the optimizer and get the response
        if options.control.case{options.control.mode}.logic.bins == 1 % binned output
            curBaseline   = Get_baseline_vals(curData,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case); % simulate the TS of baseline response
            curAll_binned = BinData(curData,options.control.case{options.control.mode}.logic.VbinSize,options.control.case{options.control.mode}.logic.PrbinSize); % bin the data according to sicretization
            Sorted_bins = Get_sorted_weighted_bins(curBaseline,curAll_binned) ; % Sort bins according to their influence on damage and revenue
            % Feed sorted bins and baseline vals to optimizer and get P per bin. It gives time series output no bins!
            [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib},options.control.case{options.control.mode}.FunSettings.wohler] = feval(options.control.case{options.control.mode}.logic.Prat,Sorted_bins,curBaseline,curAll_binned,Surrogate,curData,options.files.surrogate.case,options.control.case{options.control.mode}.logic.bins,options.control.case{options.control.mode}.FunSettings);   % Feed sorted bins and baseline vals to optimizer and get P per bin
        else % hourly values
            % Feed timeseries and baseline vals to optimizer
            curBaseline   = Get_baseline_vals(curData,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case); % simulate the TS of baseline response
            [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib},options.control.case{options.control.mode}.FunSettings.wohler] = feval(options.control.case{options.control.mode}.logic.Prat,[],curBaseline,[],Surrogate,curData,options.files.surrogate.case,options.control.case{options.control.mode}.logic.bins,options.control.case{options.control.mode}.FunSettings);   % Feed sorted bins and baseline vals to optimizer and get P per bin
        end

    elseif options.control.mode  == 6
        %    - feed the block to the optimizer and get the response
        if options.control.case{options.control.mode}.logic.bins == 1 % binned output
            % create necessary iputs for the bin optimization: discetization V-price bins/sorting according
            %  to their contribution to baseline cumulative values
            curBaseline   = Get_baseline_vals(curData,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case); % simulate the TS of baseline response
            curAll_binned = BinData(curData,options.control.case{options.control.mode}.logic.VbinSize,options.control.case{options.control.mode}.logic.PrbinSize); % bin the data according to sicretization
            Sorted_bins = Get_sorted_weighted_bins(curBaseline,curAll_binned) ; % Sort bins according to their influence on damage and revenue
            % Feed sorted bins and baseline vals to optimizer and get P per bin. It gives time series output no bins!
            [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib},options.control.case{options.control.mode}.FunSettings.wohler] = feval(options.control.case{options.control.mode}.logic.Prat,Sorted_bins,curBaseline,curAll_binned,Surrogate,curData,options.files.surrogate.case,options.control.case{options.control.mode}.logic.bins,options.control.case{options.control.mode}.FunSettings);   % Feed sorted bins and baseline vals to optimizer and get P per bin
        else % hourly values
            % Feed timeseries and baseline vals to optimizer
            curBaseline   = Get_baseline_vals(curData,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case); % simulate the TS of baseline response
            [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib},options.control.case{options.control.mode}.FunSettings.wohler] = feval(options.control.case{options.control.mode}.logic.Prat,[],curBaseline,[],Surrogate,curData,options.files.surrogate.case,options.control.case{options.control.mode}.logic.bins,options.control.case{options.control.mode}.FunSettings);   % Feed sorted bins and baseline vals to optimizer and get P per bin
            if options.control.case{options.control.mode}.FunSettings.optimizer{ib}.fval>=0
                curOpt_out(:) = options.control.case{options.control.mode}.logic.baseP ;
            end
        end

    elseif  options.control.mode  == 7
        % here I also need to feed the previous block or part of it to the next step
        curBaseline   = Get_baseline_vals(curData,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case); % simulate the TS of baseline response
        curAll_binned = BinData(curData,curBaseline,options.control.case{options.control.mode}.logic.VbinSize,options.control.case{options.control.mode}.logic.PrbinSize); % bin the data according to sicretization
        Sorted_bins = Get_sorted_weighted_bins(curBaseline,curAll_binned) ; %#ok<NASGU> % Sort bins according to their influence on damage and revenue

    elseif  options.control.mode  == 9
        curBaseline   = Get_baseline_vals_Weibull(curData,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case,options.control.case{options.control.mode}.data.Weibull,options.control.case{options.control.mode}.logic.VbinSize,options.control.case{options.control.mode}.FunSettings.ConsTI); % simulate the TS of baseline response
        % Feed values to optimizer and get P per time step and IBC/shut thresholds
        [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib},options.control.case{options.control.mode}.FunSettings.wohler] = feval(options.control.case{options.control.mode}.logic.Prat,options.control.case{options.control.mode}.logic.VbinSize,...
            options.control.case{options.control.mode}.data.Weibull,curBaseline,options.control.case{options.control.mode}.IBC,...
            Surrogate,curData,options.files.surrogate.case,options.control.case{9}.logic.baseP,options.control.case{options.control.mode}.FunSettings);

    elseif  options.control.mode  == 10
        % basicaly reads in the old result and assigns P values and IBC/shut thresholds
        [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib}] = feval(options.control.case{options.control.mode}.logic.Prat,...
            curData,options.control.case{options.control.mode}.FunSettings);


    elseif options.control.mode  == 11
        %    - feed the block to the optimizer and get the response
        if options.control.case{options.control.mode}.logic.bins == 1 % binned output
            % create necessary iputs for the bin optimization: discetization V-price bins/sorting according
            %  to their contribution to baseline cumulative values
            curBaseline   = Get_baseline_vals(curData_unc,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case); % simulate the TS of baseline response
            curAll_binned = BinData(curData_unc,options.control.case{options.control.mode}.logic.VbinSize,options.control.case{options.control.mode}.logic.PrbinSize); % bin the data according to sicretization
            Sorted_bins = Get_sorted_weighted_bins(curBaseline,curAll_binned) ; % Sort bins according to their influence on damage and revenue
            % Feed sorted bins and baseline vals to optimizer and get P per bin. It gives time series output no bins!
            [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib},options.control.case{options.control.mode}.FunSettings.wohler] = feval(options.control.case{options.control.mode}.logic.Prat,Sorted_bins,curBaseline,curAll_binned,Surrogate,curData_unc,options.files.surrogate.case,options.control.case{options.control.mode}.logic.bins,options.control.case{options.control.mode}.FunSettings);   % Feed sorted bins and baseline vals to optimizer and get P per bin
        else % hourly values
            % Feed timeseries and baseline vals to optimizer
            curBaseline   = Get_baseline_vals(curData_unc,Surrogate{1},options.control.case{options.control.mode}.logic.baseP,options.files.surrogate.case); % simulate the TS of baseline response
            [curOpt_out,options.control.case{options.control.mode}.FunSettings.optimizer{ib},options.control.case{options.control.mode}.FunSettings.wohler] = feval(options.control.case{options.control.mode}.logic.Prat,[],curBaseline,[],Surrogate,curData_unc,options.files.surrogate.case,options.control.case{options.control.mode}.logic.bins,options.control.case{options.control.mode}.FunSettings);   % Feed sorted bins and baseline vals to optimizer and get P per bin
            if options.control.case{options.control.mode}.FunSettings.optimizer{ib}.fval>=0
                curOpt_out(:) = options.control.case{options.control.mode}.logic.baseP ;
            end
        end
    end

    % looping hours of the current block and assign power level and IBC according to optimizer ot predefined logic
    for iT = 1:length(curData.Price)
        time_cnt = time_cnt+1;
        cur.price = curData.Price(iT);
        cur.V     = curData.V(iT);
        cur.TI    = curData.TI(iT);
        cur.flagoff = 0;  % signal that wind turbine is not operating

        % switch logic
        % if outside of operating range or negative prices --> shut down
        if cur.price<=0 || cur.V<4 || cur.V>24
            cur.flagoff = 1;
            cur.Sur.response =[];
            cur.cntr    = 'ShutDown';
            cur.P       = 0;
            cur.IBC     = 0;

            % 1:constant
        elseif options.control.mode == 1
            cur.P   = options.control.case{1}.Prat;
            if cur.V<options.control.IBC_generalV
                cur.IBC = 0;
            else
                cur.IBC = options.control.case{1}.IBC;
            end
            cur.cntr = 'Normal';

            % 2:shut down only
        elseif options.control.mode == 2
            cur.P   = options.control.case{2}.Prat;
            cur.off = feval(options.control.case{options.control.mode}.logic.shut,cur.V,cur.price,cur.TI,PLcurve,options.control.case{options.control.mode}.FunSettings);
            if cur.off == 1
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                cur.P   = options.control.case{options.control.mode}.Prat;
                cur.IBC = options.control.case{options.control.mode}.IBC;
                cur.cntr = 'Normal';
            end
            if cur.V<options.control.IBC_generalV
                cur.IBC = 0;
            end

            % 3:IBC only
        elseif options.control.mode == 3
            cur.P  = options.control.case{options.control.mode}.Prat;
            if cur.V<options.control.IBC_generalV
                cur.IBC = 0 ;
            else
                cur.IBC = feval(options.control.case{options.control.mode}.logic.IBC,cur.V,cur.price,cur.TI,options.control.case{options.control.mode}.FunSettings);
            end
            cur.cntr = 'Normal';

            % 4:shut down+IBC
        elseif options.control.mode == 4
            cur.off = feval(options.control.case{options.control.mode}.logic.shut,cur.V,cur.price,cur.TI,PLcurve,options.control.case{options.control.mode}.FunSettings);
            if cur.off == 1
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                if cur.V<options.control.IBC_generalV
                    cur.IBC = 0;
                else
                    cur.IBC = feval(options.control.case{options.control.mode}.logic.IBC,cur.V,cur.price,cur.TI,options.control.case{options.control.mode}.FunSettings);
                end
                cur.P   = options.control.case{options.control.mode}.Prat;
                cur.cntr = 'Normal';
            end

            % 8:shut down+IBC+Boost
        elseif options.control.mode == 8
            cur.off = feval(options.control.case{options.control.mode}.logic.shut,cur.V,cur.price,cur.TI,PLcurve,options.control.case{options.control.mode}.FunSettings);
            if cur.off == 1
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                if cur.V<options.control.IBC_generalV
                    cur.IBC = 0;
                    cur.P   = options.control.case{options.control.mode}.Prat;
                else
                    cur.P = feval(options.control.case{options.control.mode}.logic.Boost,cur.V,cur.price,cur.TI,PLcurve,options.control.case{options.control.mode}.FunSettings);
                    cur.IBC = feval(options.control.case{options.control.mode}.logic.IBC,cur.V,cur.P,cur.TI,options.control.case{options.control.mode}.FunSettings);
                end
                cur.cntr = 'Normal';
            end

            % 5:Rating with horizon
        elseif options.control.mode == 5
            cur.P = curOpt_out(iT);
            cur.cntr = 'Normal';
            if cur.V<options.control.IBC_generalV
                cur.IBC = 0;
            else
                cur.IBC = options.control.case{options.control.mode}.IBC;
            end

            % 6:Rating+shutdown+IBC with horizon
        elseif options.control.mode == 6
            % check if we are in the shut threshold. To speed up we use the precalculated Power curve
            cur.PredPow = interp2(PLcurve.Pvec,PLcurve.Vvec,PLcurve.power,curOpt_out(iT),cur.V)/1000;
            if cur.PredPow*cur.price<options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.Shut.Rev
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                cur.cntr = 'Normal';
                cur.P = curOpt_out(iT);
                if cur.V>options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.IBC.WSP
                    cur.IBC = 1;
                else
                    cur.IBC = 0;
                end
            end

            % 7: Rating+shutdown+IBC with horizon and past information
        elseif options.control.mode == 7
            cur.funIn = [cur.price cur.V cur.TI];
            cur.off = feval(options.control.case{options.control.mode}.logic.shut,cur.funIn);
            if cur.off == 1
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                cur.funIn2 = [cur.price cur.V cur.TI cur.P];
                if cur.V<options.control.IBC_generalV
                    cur.IBC = 0;
                else
                    cur.IBC = feval(options.control.case{options.control.mode}.logic.shut,cur.funIn);
                end
                if options.control.case{options.control.mode}.logic.bins == 0
                    cur.P = curOpt_out(curData.Time{1, 2}(iT));
                elseif options.control.case{options.control.mode}.logic.bins == 1
                    cur.P = curOpt_out(cur);
                end
            end

            % 9:Rating+shutdown+IBC based on wsp probability for constant prices
        elseif options.control.mode == 9
            % check if we are in the shut threshold. To speed up we use the precalculated Power curve
            cur.PredPow = interp2(PLcurve.Pvec,PLcurve.Vvec,PLcurve.power,curOpt_out(iT),cur.V)/1000;
            if cur.PredPow*cur.price<options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.Shut.Rev
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                cur.cntr = 'Normal';
                cur.P = curOpt_out(iT);
                if cur.V>options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.IBC.WSP
                    cur.IBC = 1;
                else
                    cur.IBC = 0;
                end
            end

        elseif options.control.mode == 10
            % check if we are in the shut threshold. To speed up we use the precalculated Power curve
            cur.PredPow = interp2(PLcurve.Pvec,PLcurve.Vvec,PLcurve.power,curOpt_out(iT),cur.V)/1000;
            if cur.PredPow*cur.price<options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.Shut.Rev
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                cur.cntr = 'Normal';
                cur.P = curOpt_out(iT);
                if cur.V>options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.IBC.WSP
                    cur.IBC = 1;
                else
                    cur.IBC = 0;
                end
            end

        elseif options.control.mode == 11
            % check if we are in the shut threshold. To speed up we use the precalculated Power curve
            cur.PredPow = interp2(PLcurve.Pvec,PLcurve.Vvec,PLcurve.power,curOpt_out(iT),cur.V)/1000;
            if cur.PredPow*cur.price<options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.Shut.Rev
                cur.cntr = 'ShutDown';
                cur.P    = 0;
                cur.IBC  = 0;
            else
                cur.cntr = 'Normal';
                cur.P = curOpt_out(iT);
                if cur.V>options.control.case{options.control.mode}.FunSettings.optimizer{ib}.Thresholds.IBC.WSP
                    cur.IBC = 1;
                else
                    cur.IBC = 0;
                end
            end
        end

        % Assign instantaneous values for the block
        BlockVals.num(iT,:) = [cur.P,cur.IBC];
        BlockVals.Cntr{iT,:} = cur.cntr;
        clear cur
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign instantaneous and cumulative values based on current P,wsp,TI, IBC
    BlockVals.all = [curData.V,curData.TI,BlockVals.num(:,1),BlockVals.num(:,2)]; % gather all inputs together
    %    find the time steps where the turbine operates and separate in  Normal, IBC and shut
    BlockVals.ind.NoIBC = find(BlockVals.all(:,3)~=0 & BlockVals.all(:,4)==0);
    BlockVals.ind.IBC = find(BlockVals.all(:,3)~=0 & BlockVals.all(:,4)==1);
    BlockVals.ind.Shut = find(BlockVals.all(:,3)==0);

    %Probe surrogate for all IBC and normal operation points of the block
    if options.files.surrogate.case == 1
        BlockVals.Sur.NoIBC = Get_Points_From_Surrogate_noLoad(BlockVals.all(BlockVals.ind.NoIBC,1),BlockVals.all(BlockVals.ind.NoIBC,2),BlockVals.all(BlockVals.ind.NoIBC,3),Surrogate{1,1},'spline') ;
        BlockVals.Sur.IBC = Get_Points_From_Surrogate_noLoad(BlockVals.all(BlockVals.ind.IBC,1),BlockVals.all(BlockVals.ind.IBC,2),BlockVals.all(BlockVals.ind.IBC,3),Surrogate{1,2},'spline') ;
    elseif options.files.surrogate.case  == 2
        BlockVals.Sur.NoIBC = Get_Points_From_Surrogate_GPR(BlockVals.all(BlockVals.ind.NoIBC,1),BlockVals.all(BlockVals.ind.NoIBC,2),BlockVals.all(BlockVals.ind.NoIBC,3),Surrogate{1,1});
        BlockVals.Sur.IBC = Get_Points_From_Surrogate_GPR(BlockVals.all(BlockVals.ind.IBC,1),BlockVals.all(BlockVals.ind.IBC,2),BlockVals.all(BlockVals.ind.IBC,3),Surrogate{1,2});
    end

    %Assign the values to the general output
    if ib==1
        Output = assignInstMetrics([],curData,BlockVals,ib);
        Output = assignCumMetrics(Output,ib);
    else
        Output = assignInstMetrics(Output,curData,BlockVals,ib);
        Output = assignCumMetrics(Output,ib);
    end
    clear  curData curData_unc curOpt_out curBaseline curAll_binned Sorted_bins BlockVals
end

%% Metrics and comparison to baseline
Output.metrics.CapacityFactor = Output.cum.Energy(end)/(10^4*time_cnt);
Output.metrics.IBCactivation.hours = sum(Output.inst.IBC(:)==1);
Output.metrics.IBCactivation.Perc = 100*sum(Output.inst.IBC(:)==1)/time_cnt;
Output.metrics.ShutDown.hours = sum(Output.inst.Prat(:)==0);
Output.metrics.ShutDown.Perc = sum(Output.inst.Prat(:)==0)/time_cnt;
% Output.metrics.ControlActivity = getCntPercentage(Output.inst.Cntr);

if flags.CompareToBaseline ==1
    Baseline = load(options.files.BaselineResult);
    Output.metrics.CompToBaseCum = CompareToBaselineCum(Output,Baseline,time_cnt);
end

Output.options.control = options.control.case{options.control.mode};

% if exist('FunSettings','var')
%     Output.options.control.FunSettings = options.control.case{options.control.mode}.FunSettings;
% end
Output.options.files   = options.files;
Output.options.const_price.activate = options.control.const_price.activate;
Output.options.const_price.val = options.control.const_price.val ;
Output.V     = Data_TS.V;
Output.Price = Data_TS.Price;
Output.TI    = Data_TS.TI;
Output.Time  = Data_TS.Time;

%saving
if flags.save.on == 1
    save(flags.save.savematname,'Output');
end
%%%%%%%%%%%%%%%%%% END OF MAIN CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%
% DirtyPlotting;

%% Functions
function Output = assignInstMetrics(Output,curData,BlockVals,ib)
if ib==1
    if any(BlockVals.Sur.NoIBC.Power.mean <0)
        BlockVals.Sur.NoIBC.Power.mean(BlockVals.Sur.NoIBC.Power.mean <0) = 10;
    end
    % general
    Output.inst.Cntr = BlockVals.Cntr;
    Output.inst.Prat = BlockVals.all(:,3);
    Output.inst.IBC  = BlockVals.all(:,4);

    Output.inst.Energy(BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.Power.mean ;
    Output.inst.Energy(BlockVals.ind.IBC,1)   = BlockVals.Sur.IBC.Power.mean ;
    Output.inst.Energy(BlockVals.ind.Shut,1)  = 0 ;
    Output.inst.Energy(BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.Power.std ;
    Output.inst.Energy(BlockVals.ind.IBC,2)  = BlockVals.Sur.IBC.Power.std ;
    Output.inst.Energy(BlockVals.ind.Shut,2) = 0 ;
    Output.inst.Revenue(:,1) = curData.Price.*Output.inst.Energy(:,1)/1000 ;
    Output.inst.Revenue(:,2) = curData.Price.*Output.inst.Energy(:,2)/1000 ;

    Output.inst.BlPitchTrav(BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.BlPitchTrav.mean ;
    Output.inst.BlPitchTrav(BlockVals.ind.IBC,1)   = BlockVals.Sur.IBC.BlPitchTrav.mean ;
    Output.inst.BlPitchTrav(BlockVals.ind.Shut,1)  = 0 ;
    Output.inst.BlPitchTrav(BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.BlPitchTrav.std ;
    Output.inst.BlPitchTrav(BlockVals.ind.IBC,2)   = BlockVals.Sur.IBC.Power.std ;
    Output.inst.BlPitchTrav(BlockVals.ind.Shut,2)  = 0 ;

    Output.inst.BlPitchSTD(BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.BlPitchSTD.mean ;
    Output.inst.BlPitchSTD(BlockVals.ind.IBC,1)   = BlockVals.Sur.IBC.BlPitchSTD.mean ;
    Output.inst.BlPitchSTD(BlockVals.ind.Shut,1)  = 0 ;
    Output.inst.BlPitchSTD(BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.BlPitchSTD.std ;
    Output.inst.BlPitchSTD(BlockVals.ind.IBC,2)   = BlockVals.Sur.IBC.BlPitchSTD.std ;
    Output.inst.BlPitchSTD(BlockVals.ind.Shut,2)  = 0 ;

    Output.inst.GenTqSTD(BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.GenTqSTD.mean ;
    Output.inst.GenTqSTD(BlockVals.ind.IBC,1)   = BlockVals.Sur.IBC.GenTqSTD.mean ;
    Output.inst.GenTqSTD(BlockVals.ind.Shut,1)  = 0 ;
    Output.inst.GenTqSTD(BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.GenTqSTD.std ;
    Output.inst.GenTqSTD(BlockVals.ind.IBC,2)   = BlockVals.Sur.IBC.GenTqSTD.std ;
    Output.inst.GenTqSTD(BlockVals.ind.Shut,2)  = 0 ;

    Output.inst.GenSpeedSTD(BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.GenSpeedSTD.mean ;
    Output.inst.GenSpeedSTD(BlockVals.ind.IBC,1)   = BlockVals.Sur.IBC.GenSpeedSTD.mean ;
    Output.inst.GenSpeedSTD(BlockVals.ind.Shut,1)  = 0 ;
    Output.inst.GenSpeedSTD(BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.GenSpeedSTD.std ;
    Output.inst.GenSpeedSTD(BlockVals.ind.IBC,2)  = BlockVals.Sur.IBC.GenSpeedSTD.std ;
    Output.inst.GenSpeedSTD(BlockVals.ind.Shut,2) = 0 ;

    %DEL
    Output.inst.DEL.TBMx(BlockVals.ind.NoIBC,:)  = [BlockVals.Sur.NoIBC.TBMx.mean,BlockVals.Sur.NoIBC.TBMx.std] ;
    Output.inst.DEL.TBMx(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TBMx.mean,BlockVals.Sur.IBC.TBMx.std] ;
    Output.inst.DEL.TBMx(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TBMy(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TBMy.mean,BlockVals.Sur.NoIBC.TBMy.std] ;
    Output.inst.DEL.TBMy(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TBMy.mean,BlockVals.Sur.IBC.TBMy.std] ;
    Output.inst.DEL.TBMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TBMz(BlockVals.ind.NoIBC,:)  = [BlockVals.Sur.NoIBC.TBMz.mean,BlockVals.Sur.NoIBC.TBMz.std] ;
    Output.inst.DEL.TBMz(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TBMz.mean,BlockVals.Sur.IBC.TBMz.std] ;
    Output.inst.DEL.TBMz(BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DEL.BRMx(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRMx.mean,BlockVals.Sur.NoIBC.BRMx.std] ;
    Output.inst.DEL.BRMx(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRMx.mean,BlockVals.Sur.IBC.BRMx.std] ;
    Output.inst.DEL.BRMx(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.BRMy(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRMy.mean,BlockVals.Sur.NoIBC.BRMy.std] ;
    Output.inst.DEL.BRMy(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRMy.mean,BlockVals.Sur.IBC.BRMy.std] ;
    Output.inst.DEL.BRMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.BRMz(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRMz.mean,BlockVals.Sur.NoIBC.BRMz.std] ;
    Output.inst.DEL.BRMz(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRMz.mean,BlockVals.Sur.IBC.BRMz.std] ;
    Output.inst.DEL.BRMz(BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DEL.BROop(BlockVals.ind.NoIBC,:)= [BlockVals.Sur.NoIBC.BROop.mean,BlockVals.Sur.NoIBC.BROop.std] ;
    Output.inst.DEL.BROop(BlockVals.ind.IBC,:)  = [BlockVals.Sur.IBC.BROop.mean,BlockVals.Sur.IBC.BROop.std] ;
    Output.inst.DEL.BROop(BlockVals.ind.Shut,:) = 0  ;
    Output.inst.DEL.BRIp(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRIp.mean,BlockVals.Sur.NoIBC.BRIp.std] ;
    Output.inst.DEL.BRIp(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRIp.mean,BlockVals.Sur.IBC.BRIp.std] ;
    Output.inst.DEL.BRIp(BlockVals.ind.Shut,:)  = 0  ;


    Output.inst.DEL.TTMx(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TTMx.mean,BlockVals.Sur.NoIBC.TTMx.std] ;
    Output.inst.DEL.TTMx(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TTMx.mean,BlockVals.Sur.IBC.TTMx.std] ;
    Output.inst.DEL.TTMx(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TTMy(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TTMy.mean,BlockVals.Sur.NoIBC.TTMy.std] ;
    Output.inst.DEL.TTMy(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TTMy.mean,BlockVals.Sur.IBC.TTMy.std] ;
    Output.inst.DEL.TTMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TTMz(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TTMz.mean,BlockVals.Sur.NoIBC.TTMz.std] ;
    Output.inst.DEL.TTMz(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TTMz.mean,BlockVals.Sur.IBC.TTMz.std] ;
    Output.inst.DEL.TTMz(BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DEL.LSSMy(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.LSSMy.mean,BlockVals.Sur.NoIBC.LSSMy.std] ;
    Output.inst.DEL.LSSMy(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.LSSMy.mean,BlockVals.Sur.IBC.LSSMy.std] ;
    Output.inst.DEL.LSSMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.LSSMz(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.LSSMz.mean,BlockVals.Sur.NoIBC.LSSMz.std] ;
    Output.inst.DEL.LSSMz(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.LSSMz.mean,BlockVals.Sur.IBC.LSSMz.std] ;
    Output.inst.DEL.LSSMz(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.LSSTq(BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.LSSTq.mean,BlockVals.Sur.NoIBC.LSSTq.std] ;
    Output.inst.DEL.LSSTq(BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.LSSTq.mean,BlockVals.Sur.IBC.LSSTq.std] ;
    Output.inst.DEL.LSSTq(BlockVals.ind.Shut,:)  = 0  ;
    % Damage
    Output.inst.DAM.TBMx(BlockVals.ind.NoIBC,:)  = [(BlockVals.Sur.NoIBC.TBMx.mean.^4)*3600,(BlockVals.Sur.NoIBC.TBMx.std.^4)*3600] ;
    Output.inst.DAM.TBMx(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TBMx.mean.^4)*3600,(BlockVals.Sur.IBC.TBMx.std.^4)*3600] ;
    Output.inst.DAM.TBMx(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TBMy(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TBMy.mean.^4)*3600,(BlockVals.Sur.NoIBC.TBMy.std.^4)*3600] ;
    Output.inst.DAM.TBMy(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TBMy.mean.^4)*3600,(BlockVals.Sur.IBC.TBMy.std.^4)*3600] ;
    Output.inst.DAM.TBMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TBMz(BlockVals.ind.NoIBC,:)  = [(BlockVals.Sur.NoIBC.TBMz.mean.^4)*3600,(BlockVals.Sur.NoIBC.TBMz.std.^4)*3600] ;
    Output.inst.DAM.TBMz(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TBMz.mean.^4)*3600,(BlockVals.Sur.IBC.TBMz.std.^4)*3600] ;
    Output.inst.DAM.TBMz(BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DAM.BRMx(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRMx.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRMx.std.^10)*3600] ;
    Output.inst.DAM.BRMx(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRMx.mean.^10)*3600,(BlockVals.Sur.IBC.BRMx.std.^10)*3600] ;
    Output.inst.DAM.BRMx(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.BRMy(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRMy.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRMy.std.^10)*3600] ;
    Output.inst.DAM.BRMy(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRMy.mean.^10)*3600,(BlockVals.Sur.IBC.BRMy.std.^10)*3600] ;
    Output.inst.DAM.BRMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.BRMz(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRMz.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRMz.std.^10)*3600] ;
    Output.inst.DAM.BRMz(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRMz.mean.^10)*3600,(BlockVals.Sur.IBC.BRMz.std.^10)*3600] ;
    Output.inst.DAM.BRMz(BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DAM.BROop(BlockVals.ind.NoIBC,:)= [(BlockVals.Sur.NoIBC.BROop.mean.^10)*3600,(BlockVals.Sur.NoIBC.BROop.std.^10)*3600] ;
    Output.inst.DAM.BROop(BlockVals.ind.IBC,:)  = [(BlockVals.Sur.IBC.BROop.mean.^10)*3600,(BlockVals.Sur.IBC.BROop.std.^10)*3600] ;
    Output.inst.DAM.BROop(BlockVals.ind.Shut,:) = 0  ;
    Output.inst.DAM.BRIp(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRIp.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRIp.std.^10)*3600] ;
    Output.inst.DAM.BRIp(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRIp.mean.^10)*3600,(BlockVals.Sur.IBC.BRIp.std.^10)*3600] ;
    Output.inst.DAM.BRIp(BlockVals.ind.Shut,:)  = 0  ;


    Output.inst.DAM.TTMx(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TTMx.mean.^4)*3600,(BlockVals.Sur.NoIBC.TTMx.std.^4)*3600] ;
    Output.inst.DAM.TTMx(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TTMx.mean.^4)*3600,(BlockVals.Sur.IBC.TTMx.std.^4)*3600] ;
    Output.inst.DAM.TTMx(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TTMy(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TTMy.mean.^4)*3600,(BlockVals.Sur.NoIBC.TTMy.std.^4)*3600] ;
    Output.inst.DAM.TTMy(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TTMy.mean.^4)*3600,(BlockVals.Sur.IBC.TTMy.std.^4)*3600] ;
    Output.inst.DAM.TTMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TTMz(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TTMz.mean.^4)*3600,(BlockVals.Sur.NoIBC.TTMz.std.^4)*3600] ;
    Output.inst.DAM.TTMz(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TTMz.mean.^4)*3600,(BlockVals.Sur.IBC.TTMz.std.^4)*3600] ;
    Output.inst.DAM.TTMz(BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DAM.LSSMy(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.LSSMy.mean.^4)*3600,(BlockVals.Sur.NoIBC.LSSMy.std.^4)*3600] ;
    Output.inst.DAM.LSSMy(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.LSSMy.mean.^4)*3600,(BlockVals.Sur.IBC.LSSMy.std.^4)*3600] ;
    Output.inst.DAM.LSSMy(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.LSSMz(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.LSSMz.mean.^4)*3600,(BlockVals.Sur.NoIBC.LSSMz.std.^4)*3600] ;
    Output.inst.DAM.LSSMz(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.LSSMz.mean.^4)*3600,(BlockVals.Sur.IBC.LSSMz.std.^4)*3600] ;
    Output.inst.DAM.LSSMz(BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.LSSTq(BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.LSSTq.mean.^4)*3600,(BlockVals.Sur.NoIBC.LSSTq.std.^4)*3600] ;
    Output.inst.DAM.LSSTq(BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.LSSTq.mean.^4)*3600,(BlockVals.Sur.IBC.LSSTq.std.^4)*3600] ;
    Output.inst.DAM.LSSTq(BlockVals.ind.Shut,:)  = 0  ;


else
    if any(BlockVals.Sur.NoIBC.Power.mean <0)
        BlockVals.Sur.NoIBC.Power.mean(BlockVals.Sur.NoIBC.Power.mean <0) = 10;
    end
    PrevL = length(Output.inst.Cntr);
    %general
    Output.inst.Cntr = [Output.inst.Cntr; BlockVals.Cntr];
    Output.inst.Prat = [Output.inst.Prat; BlockVals.all(:,3)];
    Output.inst.IBC  = [Output.inst.IBC; BlockVals.all(:,4)];
    Output.inst.Energy(PrevL+BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.Power.mean ;
    Output.inst.Energy(PrevL+BlockVals.ind.IBC,1) = BlockVals.Sur.IBC.Power.mean ;
    Output.inst.Energy(PrevL+BlockVals.ind.Shut,1) = 0 ;
    Output.inst.Energy(PrevL+BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.Power.std ;
    Output.inst.Energy(PrevL+BlockVals.ind.IBC,2) = BlockVals.Sur.IBC.Power.std ;
    Output.inst.Energy(PrevL+BlockVals.ind.Shut,2) = 0 ;
    Output.inst.Revenue = [Output.inst.Revenue;curData.Price.*Output.inst.Energy(PrevL+1:end,1)/1000,curData.Price.*Output.inst.Energy(PrevL+1:end,2)/1000] ;


    Output.inst.BlPitchTrav(PrevL+BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.BlPitchTrav.mean ;
    Output.inst.BlPitchTrav(PrevL+BlockVals.ind.IBC,1) = BlockVals.Sur.IBC.BlPitchTrav.mean ;
    Output.inst.BlPitchTrav(PrevL+BlockVals.ind.Shut,1) = 0 ;
    Output.inst.BlPitchTrav(PrevL+BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.BlPitchTrav.std ;
    Output.inst.BlPitchTrav(PrevL+BlockVals.ind.IBC,2) = BlockVals.Sur.IBC.Power.std ;
    Output.inst.BlPitchTrav(PrevL+BlockVals.ind.Shut,2) = 0 ;

    Output.inst.BlPitchSTD(PrevL+BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.BlPitchSTD.mean ;
    Output.inst.BlPitchSTD(PrevL+BlockVals.ind.IBC,1) = BlockVals.Sur.IBC.BlPitchSTD.mean ;
    Output.inst.BlPitchSTD(PrevL+BlockVals.ind.Shut,1) = 0 ;
    Output.inst.BlPitchSTD(PrevL+BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.BlPitchSTD.std ;
    Output.inst.BlPitchSTD(PrevL+BlockVals.ind.IBC,2) = BlockVals.Sur.IBC.BlPitchSTD.std ;
    Output.inst.BlPitchSTD(PrevL+BlockVals.ind.Shut,2) = 0 ;

    Output.inst.GenTqSTD(PrevL+BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.GenTqSTD.mean ;
    Output.inst.GenTqSTD(PrevL+BlockVals.ind.IBC,1) = BlockVals.Sur.IBC.GenTqSTD.mean ;
    Output.inst.GenTqSTD(PrevL+BlockVals.ind.Shut,1) = 0 ;
    Output.inst.GenTqSTD(PrevL+BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.GenTqSTD.std ;
    Output.inst.GenTqSTD(PrevL+BlockVals.ind.IBC,2) = BlockVals.Sur.IBC.GenTqSTD.std ;
    Output.inst.GenTqSTD(PrevL+BlockVals.ind.Shut,2) = 0 ;

    Output.inst.GenSpeedSTD(PrevL+BlockVals.ind.NoIBC,1) = BlockVals.Sur.NoIBC.GenSpeedSTD.mean ;
    Output.inst.GenSpeedSTD(PrevL+BlockVals.ind.IBC,1) = BlockVals.Sur.IBC.GenSpeedSTD.mean ;
    Output.inst.GenSpeedSTD(PrevL+BlockVals.ind.Shut,1) = 0 ;
    Output.inst.GenSpeedSTD(PrevL+BlockVals.ind.NoIBC,2) = BlockVals.Sur.NoIBC.GenSpeedSTD.std ;
    Output.inst.GenSpeedSTD(PrevL+BlockVals.ind.IBC,2) = BlockVals.Sur.IBC.GenSpeedSTD.std ;
    Output.inst.GenSpeedSTD(PrevL+BlockVals.ind.Shut,2) = 0 ;
    % DEL
    Output.inst.DEL.TBMx(PrevL+BlockVals.ind.NoIBC,:)  = [BlockVals.Sur.NoIBC.TBMx.mean,BlockVals.Sur.NoIBC.TBMx.std] ;
    Output.inst.DEL.TBMx(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TBMx.mean,BlockVals.Sur.IBC.TBMx.std] ;
    Output.inst.DEL.TBMx(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TBMy(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TBMy.mean,BlockVals.Sur.NoIBC.TBMy.std] ;
    Output.inst.DEL.TBMy(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TBMy.mean,BlockVals.Sur.IBC.TBMy.std] ;
    Output.inst.DEL.TBMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TBMz(PrevL+BlockVals.ind.NoIBC,:)  = [BlockVals.Sur.NoIBC.TBMz.mean,BlockVals.Sur.NoIBC.TBMz.std] ;
    Output.inst.DEL.TBMz(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TBMz.mean,BlockVals.Sur.IBC.TBMz.std] ;
    Output.inst.DEL.TBMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DEL.BRMx(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRMx.mean,BlockVals.Sur.NoIBC.BRMx.std] ;
    Output.inst.DEL.BRMx(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRMx.mean,BlockVals.Sur.IBC.BRMx.std] ;
    Output.inst.DEL.BRMx(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.BRMy(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRMy.mean,BlockVals.Sur.NoIBC.BRMy.std] ;
    Output.inst.DEL.BRMy(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRMy.mean,BlockVals.Sur.IBC.BRMy.std] ;
    Output.inst.DEL.BRMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.BRMz(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRMz.mean,BlockVals.Sur.NoIBC.BRMz.std] ;
    Output.inst.DEL.BRMz(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRMz.mean,BlockVals.Sur.IBC.BRMz.std] ;
    Output.inst.DEL.BRMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DEL.BROop(PrevL+BlockVals.ind.NoIBC,:)= [BlockVals.Sur.NoIBC.BROop.mean,BlockVals.Sur.NoIBC.BROop.std] ;
    Output.inst.DEL.BROop(PrevL+BlockVals.ind.IBC,:)  = [BlockVals.Sur.IBC.BROop.mean,BlockVals.Sur.IBC.BROop.std] ;
    Output.inst.DEL.BROop(PrevL+BlockVals.ind.Shut,:) = 0  ;
    Output.inst.DEL.BRIp(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.BRIp.mean,BlockVals.Sur.NoIBC.BRIp.std] ;
    Output.inst.DEL.BRIp(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.BRIp.mean,BlockVals.Sur.IBC.BRIp.std] ;
    Output.inst.DEL.BRIp(PrevL+BlockVals.ind.Shut,:)  = 0  ;


    Output.inst.DEL.TTMx(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TTMx.mean,BlockVals.Sur.NoIBC.TTMx.std] ;
    Output.inst.DEL.TTMx(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TTMx.mean,BlockVals.Sur.IBC.TTMx.std] ;
    Output.inst.DEL.TTMx(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TTMy(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TTMy.mean,BlockVals.Sur.NoIBC.TTMy.std] ;
    Output.inst.DEL.TTMy(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TTMy.mean,BlockVals.Sur.IBC.TTMy.std] ;
    Output.inst.DEL.TTMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.TTMz(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.TTMz.mean,BlockVals.Sur.NoIBC.TTMz.std] ;
    Output.inst.DEL.TTMz(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.TTMz.mean,BlockVals.Sur.IBC.TTMz.std] ;
    Output.inst.DEL.TTMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DEL.LSSMy(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.LSSMy.mean,BlockVals.Sur.NoIBC.LSSMy.std] ;
    Output.inst.DEL.LSSMy(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.LSSMy.mean,BlockVals.Sur.IBC.LSSMy.std] ;
    Output.inst.DEL.LSSMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.LSSMz(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.LSSMz.mean,BlockVals.Sur.NoIBC.LSSMz.std] ;
    Output.inst.DEL.LSSMz(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.LSSMz.mean,BlockVals.Sur.IBC.LSSMz.std] ;
    Output.inst.DEL.LSSMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DEL.LSSTq(PrevL+BlockVals.ind.NoIBC,:) = [BlockVals.Sur.NoIBC.LSSTq.mean,BlockVals.Sur.NoIBC.LSSTq.std] ;
    Output.inst.DEL.LSSTq(PrevL+BlockVals.ind.IBC,:)   = [BlockVals.Sur.IBC.LSSTq.mean,BlockVals.Sur.IBC.LSSTq.std] ;
    Output.inst.DEL.LSSTq(PrevL+BlockVals.ind.Shut,:)  = 0  ;

    %%% Damage
    Output.inst.DAM.TBMx(PrevL+BlockVals.ind.NoIBC,:)  = [(BlockVals.Sur.NoIBC.TBMx.mean.^4)*3600,(BlockVals.Sur.NoIBC.TBMx.std.^4)*3600] ;
    Output.inst.DAM.TBMx(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TBMx.mean.^4)*3600,(BlockVals.Sur.IBC.TBMx.std.^4)*3600] ;
    Output.inst.DAM.TBMx(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TBMy(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TBMy.mean.^4)*3600,(BlockVals.Sur.NoIBC.TBMy.std.^4)*3600] ;
    Output.inst.DAM.TBMy(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TBMy.mean.^4)*3600,(BlockVals.Sur.IBC.TBMy.std.^4)*3600] ;
    Output.inst.DAM.TBMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TBMz(PrevL+BlockVals.ind.NoIBC,:)  = [(BlockVals.Sur.NoIBC.TBMz.mean.^4)*3600,(BlockVals.Sur.NoIBC.TBMz.std.^4)*3600] ;
    Output.inst.DAM.TBMz(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TBMz.mean.^4)*3600,(BlockVals.Sur.IBC.TBMz.std.^4)*3600] ;
    Output.inst.DAM.TBMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DAM.BRMx(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRMx.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRMx.std.^10)*3600] ;
    Output.inst.DAM.BRMx(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRMx.mean.^10)*3600,(BlockVals.Sur.IBC.BRMx.std.^10)*3600] ;
    Output.inst.DAM.BRMx(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.BRMy(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRMy.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRMy.std.^10)*3600] ;
    Output.inst.DAM.BRMy(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRMy.mean.^10)*3600,(BlockVals.Sur.IBC.BRMy.std.^10)*3600] ;
    Output.inst.DAM.BRMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.BRMz(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRMz.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRMz.std.^10)*3600] ;
    Output.inst.DAM.BRMz(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRMz.mean.^10)*3600,(BlockVals.Sur.IBC.BRMz.std.^10)*3600] ;
    Output.inst.DAM.BRMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DAM.BROop(PrevL+BlockVals.ind.NoIBC,:)= [(BlockVals.Sur.NoIBC.BROop.mean.^10)*3600,(BlockVals.Sur.NoIBC.BROop.std.^10)*3600] ;
    Output.inst.DAM.BROop(PrevL+BlockVals.ind.IBC,:)  = [(BlockVals.Sur.IBC.BROop.mean.^10)*3600,(BlockVals.Sur.IBC.BROop.std.^10)*3600] ;
    Output.inst.DAM.BROop(PrevL+BlockVals.ind.Shut,:) = 0  ;
    Output.inst.DAM.BRIp(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.BRIp.mean.^10)*3600,(BlockVals.Sur.NoIBC.BRIp.std.^10)*3600] ;
    Output.inst.DAM.BRIp(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.BRIp.mean.^10)*3600,(BlockVals.Sur.IBC.BRIp.std.^10)*3600] ;
    Output.inst.DAM.BRIp(PrevL+BlockVals.ind.Shut,:)  = 0  ;


    Output.inst.DAM.TTMx(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TTMx.mean.^4)*3600,(BlockVals.Sur.NoIBC.TTMx.std.^4)*3600] ;
    Output.inst.DAM.TTMx(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TTMx.mean.^4)*3600,(BlockVals.Sur.IBC.TTMx.std.^4)*3600] ;
    Output.inst.DAM.TTMx(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TTMy(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TTMy.mean.^4)*3600,(BlockVals.Sur.NoIBC.TTMy.std.^4)*3600] ;
    Output.inst.DAM.TTMy(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TTMy.mean.^4)*3600,(BlockVals.Sur.IBC.TTMy.std.^4)*3600] ;
    Output.inst.DAM.TTMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.TTMz(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.TTMz.mean.^4)*3600,(BlockVals.Sur.NoIBC.TTMz.std.^4)*3600] ;
    Output.inst.DAM.TTMz(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.TTMz.mean.^4)*3600,(BlockVals.Sur.IBC.TTMz.std.^4)*3600] ;
    Output.inst.DAM.TTMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;

    Output.inst.DAM.LSSMy(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.LSSMy.mean.^4)*3600,(BlockVals.Sur.NoIBC.LSSMy.std.^4)*3600] ;
    Output.inst.DAM.LSSMy(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.LSSMy.mean.^4)*3600,(BlockVals.Sur.IBC.LSSMy.std.^4)*3600] ;
    Output.inst.DAM.LSSMy(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.LSSMz(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.LSSMz.mean.^4)*3600,(BlockVals.Sur.NoIBC.LSSMz.std.^4)*3600] ;
    Output.inst.DAM.LSSMz(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.LSSMz.mean.^4)*3600,(BlockVals.Sur.IBC.LSSMz.std.^4)*3600] ;
    Output.inst.DAM.LSSMz(PrevL+BlockVals.ind.Shut,:)  = 0  ;
    Output.inst.DAM.LSSTq(PrevL+BlockVals.ind.NoIBC,:) = [(BlockVals.Sur.NoIBC.LSSTq.mean.^4)*3600,(BlockVals.Sur.NoIBC.LSSTq.std.^4)*3600] ;
    Output.inst.DAM.LSSTq(PrevL+BlockVals.ind.IBC,:)   = [(BlockVals.Sur.IBC.LSSTq.mean.^4)*3600,(BlockVals.Sur.IBC.LSSTq.std.^4)*3600] ;
    Output.inst.DAM.LSSTq(PrevL+BlockVals.ind.Shut,:)  = 0  ;


    clear prevL
end
end

function Output = assignCumMetrics(Output,ib)

if ib==1
    Output.cum.Energy = cumsum(Output.inst.Energy(:,1)) ;
    Output.cum.Revenue = cumsum(Output.inst.Revenue(:,1)) ;
    Output.cum.BlPitchTrav = cumsum(Output.inst.BlPitchTrav(:,1)) ;

    Output.cum.DAM.TBMx = cumsum(Output.inst.DAM.TBMx(:,1));
    Output.cum.DAM.TBMy = cumsum(Output.inst.DAM.TBMy(:,1));
    Output.cum.DAM.TBMz = cumsum(Output.inst.DAM.TBMz(:,1));

    Output.cum.DAM.TTMx = cumsum(Output.inst.DAM.TTMx(:,1));
    Output.cum.DAM.TTMy = cumsum(Output.inst.DAM.TTMy(:,1));
    Output.cum.DAM.TTMz = cumsum(Output.inst.DAM.TTMz(:,1));

    Output.cum.DAM.BRMx = cumsum(Output.inst.DAM.BRMx(:,1));
    Output.cum.DAM.BRMy = cumsum(Output.inst.DAM.BRMy(:,1));
    Output.cum.DAM.BRMz = cumsum(Output.inst.DAM.BRMz(:,1));

    Output.cum.DAM.BROop = cumsum(Output.inst.DAM.BROop(:,1));
    Output.cum.DAM.BRIp = cumsum(Output.inst.DAM.BRIp(:,1));

    Output.cum.DAM.LSSMy = cumsum(Output.inst.DAM.LSSMy(:,1));
    Output.cum.DAM.LSSMz = cumsum(Output.inst.DAM.LSSMz(:,1));
    Output.cum.DAM.LSSTq = cumsum(Output.inst.DAM.LSSTq(:,1));
else
    %     PrevL = length(Output.cum.Energy);
    Output.cum.Energy = cumsum(Output.inst.Energy(:,1));%[Output.cum.Energy(1:end-1,1);cumsum([Output.cum.Energy(end,1);Output.inst.Energy(PrevL+1:end,1)])] ;
    Output.cum.Revenue = cumsum(Output.inst.Revenue(:,1)); %[Output.cum.Revenue(1:end-1,1) ;cumsum([Output.cum.Revenue(end,1) ;Output.inst.Revenue(PrevL+1:end,1)])] ;
    Output.cum.BlPitchTrav = cumsum(Output.inst.BlPitchTrav(:,1)) ;

    Output.cum.DAM.TBMx = cumsum(Output.inst.DAM.TBMx(:,1));
    Output.cum.DAM.TBMy = cumsum(Output.inst.DAM.TBMy(:,1));
    Output.cum.DAM.TBMz = cumsum(Output.inst.DAM.TBMz(:,1));

    Output.cum.DAM.TTMx = cumsum(Output.inst.DAM.TTMx(:,1));
    Output.cum.DAM.TTMy = cumsum(Output.inst.DAM.TTMy(:,1));
    Output.cum.DAM.TTMz = cumsum(Output.inst.DAM.TTMz(:,1));

    Output.cum.DAM.BRMx = cumsum(Output.inst.DAM.BRMx(:,1));
    Output.cum.DAM.BRMy = cumsum(Output.inst.DAM.BRMy(:,1));
    Output.cum.DAM.BRMz = cumsum(Output.inst.DAM.BRMz(:,1));

    Output.cum.DAM.BROop = cumsum(Output.inst.DAM.BROop(:,1));
    Output.cum.DAM.BRIp = cumsum(Output.inst.DAM.BRIp(:,1));

    Output.cum.DAM.LSSMy = cumsum(Output.inst.DAM.LSSMy(:,1));
    Output.cum.DAM.LSSMz = cumsum(Output.inst.DAM.LSSMz(:,1));
    Output.cum.DAM.LSSTq = cumsum(Output.inst.DAM.LSSTq(:,1));

end


end

function CompToBase = CompareToBaselineCum(Output,Baseline,time_cnt)
CompToBase.Energy  = 100*(Output.cum.Energy(end)/Baseline.Output.cum.Energy(end)-1);
CompToBase.Revenue = 100*(Output.cum.Revenue(end)/Baseline.Output.cum.Revenue(end)-1);
CompToBase.TBMx = 100*(Output.cum.DAM.TBMx(end)/Baseline.Output.cum.DAM.TBMx(end)-1);
CompToBase.TBMy = 100*(Output.cum.DAM.TBMy(end)/Baseline.Output.cum.DAM.TBMy(end)-1);
CompToBase.TBMz = 100*(Output.cum.DAM.TBMz(end)/Baseline.Output.cum.DAM.TBMz(end)-1);
CompToBase.BRMx = 100*(Output.cum.DAM.BRMx(end)/Baseline.Output.cum.DAM.BRMx(end)-1);
CompToBase.BRMy = 100*(Output.cum.DAM.BRMy(end)/Baseline.Output.cum.DAM.BRMy(end)-1);
CompToBase.BRMz = 100*(Output.cum.DAM.BRMz(end)/Baseline.Output.cum.DAM.BRMz(end)-1);
CompToBase.BROop = 100*(Output.cum.DAM.BROop(end)/Baseline.Output.cum.DAM.BROop(end)-1);
CompToBase.BRIp = 100*(Output.cum.DAM.BRIp(end)/Baseline.Output.cum.DAM.BRIp(end)-1);
CompToBase.TTMx = 100*(Output.cum.DAM.TTMx(end)/Baseline.Output.cum.DAM.TTMx(end)-1);
CompToBase.TTMy = 100*(Output.cum.DAM.TTMy(end)/Baseline.Output.cum.DAM.TTMy(end)-1);
CompToBase.TTMz = 100*(Output.cum.DAM.TTMz(end)/Baseline.Output.cum.DAM.TTMz(end)-1);
CompToBase.LSSMy = 100*(Output.cum.DAM.LSSMy(end)/Baseline.Output.cum.DAM.LSSMy(end)-1);
CompToBase.LSSMz = 100*(Output.cum.DAM.LSSMz(end)/Baseline.Output.cum.DAM.LSSMz(end)-1);
CompToBase.LSSTq = 100*(Output.cum.DAM.LSSTq(end)/Baseline.Output.cum.DAM.LSSTq(end)-1);
CompToBase.BlPitchTrav = 100*(Output.cum.BlPitchTrav(end,1)/Baseline.Output.cum.BlPitchTrav(end,1)-1);
CompToBase.GenTqSTD    = 100*(mean(Output.inst.GenTqSTD(:,1))/mean(Baseline.Output.inst.GenTqSTD(:,1))-1);
CompToBase.GenSpeedSTD = 100*(mean(Output.inst.GenSpeedSTD(:,1))/mean(Baseline.Output.inst.GenSpeedSTD(:,1))-1);
CompToBase.BlPitchSTD  = 100*(mean(Output.inst.BlPitchSTD(:,1))/mean(Baseline.Output.inst.BlPitchSTD(:,1))-1);
CompToBase.ShutRelPerc  = 100*(Output.metrics.ShutDown.Perc/Baseline.Output.metrics.ShutDown.Perc-1);
CompToBase.ShutHoursDiff = Output.metrics.ShutDown.hours-Baseline.Output.metrics.ShutDown.hours;
CompToBase.ShutAbsHoursPerc = 100*CompToBase.ShutHoursDiff/time_cnt;
end

% function DELinst = assignInstDEL(cur_surResponse,instVec,flagoff)
%
% if flagoff ==1
%     cur_surResponse.TBMx.mean = 0;
%     cur_surResponse.TBMy.mean = 0;
%     cur_surResponse.TBMz.mean = 0;
%     cur_surResponse.BRMx.mean = 0;
%     cur_surResponse.BRMy.mean = 0;
%     cur_surResponse.BRMz.mean = 0;
%     cur_surResponse.BROop.mean = 0;
%     cur_surResponse.BRIp.mean = 0;
%     cur_surResponse.TTMx.mean = 0;
%     cur_surResponse.TTMy.mean = 0;
%     cur_surResponse.TTMz.mean = 0;
%     cur_surResponse.TTMz.mean = 0;
%     cur_surResponse.LSSMy.mean = 0;
%     cur_surResponse.LSSMz.mean = 0;
%     cur_surResponse.LSSTq.mean = 0;
%
%     cur_surResponse.TBMx.std = 0;
%     cur_surResponse.TBMy.std = 0;
%     cur_surResponse.TBMz.std = 0;
%     cur_surResponse.BRMx.std = 0;
%     cur_surResponse.BRMy.std = 0;
%     cur_surResponse.BRMz.std = 0;
%     cur_surResponse.BROop.std = 0;
%     cur_surResponse.BRIp.std = 0;
%     cur_surResponse.TTMx.std = 0;
%     cur_surResponse.TTMy.std = 0;
%     cur_surResponse.TTMz.std = 0;
%     cur_surResponse.TTMz.std = 0;
%     cur_surResponse.LSSMy.std = 0;
%     cur_surResponse.LSSMz.std = 0;
%     cur_surResponse.LSSTq.std = 0;
% end
%
% if isempty(instVec)
%     DELinst.TBMx.mean  = cur_surResponse.TBMx.mean ;
%     DELinst.TBMy.mean  = cur_surResponse.TBMy.mean ;
%     DELinst.TBMz.mean  = cur_surResponse.TBMz.mean ;
%     DELinst.BRMx.mean  = cur_surResponse.BRMx.mean ;
%     DELinst.BRMy.mean  = cur_surResponse.BRMy.mean ;
%     DELinst.BRMz.mean  = cur_surResponse.BRMz.mean ;
%     DELinst.BROop.mean = cur_surResponse.BROop.mean ;
%     DELinst.BRIp.mean  = cur_surResponse.BRIp.mean ;
%     DELinst.TTMx.mean  = cur_surResponse.TTMx.mean ;
%     DELinst.TTMy.mean  = cur_surResponse.TTMy.mean ;
%     DELinst.TTMz.mean  = cur_surResponse.TTMz.mean ;
%     DELinst.LSSMy.mean = cur_surResponse.LSSMy.mean ;
%     DELinst.LSSMz.mean = cur_surResponse.LSSMz.mean ;
%     DELinst.LSSTq.mean = cur_surResponse.LSSTq.mean ;
%
%     DELinst.TBMx.std  = cur_surResponse.TBMx.std ;
%     DELinst.TBMy.std  = cur_surResponse.TBMy.std ;
%     DELinst.TBMz.std  = cur_surResponse.TBMz.std ;
%     DELinst.BRMx.std  = cur_surResponse.BRMx.std ;
%     DELinst.BRMy.std  = cur_surResponse.BRMy.std ;
%     DELinst.BRMz.std  = cur_surResponse.BRMz.std ;
%     DELinst.BROop.std = cur_surResponse.BROop.std ;
%     DELinst.BRIp.std  = cur_surResponse.BRIp.std ;
%     DELinst.TTMx.std  = cur_surResponse.TTMx.std ;
%     DELinst.TTMy.std  = cur_surResponse.TTMy.std ;
%     DELinst.TTMz.std  = cur_surResponse.TTMz.std ;
%     DELinst.LSSMy.std = cur_surResponse.LSSMy.std ;
%     DELinst.LSSMz.std = cur_surResponse.LSSMz.std ;
%     DELinst.LSSTq.std = cur_surResponse.LSSTq.std ;
% else
%     DELinst.TBMx.mean  = [instVec.TBMx.mean ;cur_surResponse.TBMx.mean] ;
%     DELinst.TBMy.mean  = [instVec.TBMy.mean ;cur_surResponse.TBMy.mean] ;
%     DELinst.TBMz.mean  = [instVec.TBMz.mean ;cur_surResponse.TBMz.mean] ;
%     DELinst.BRMx.mean  = [instVec.BRMx.mean ;cur_surResponse.BRMx.mean] ;
%     DELinst.BRMy.mean  = [instVec.BRMy.mean ;cur_surResponse.BRMy.mean] ;
%     DELinst.BRMz.mean  = [instVec.BRMz.mean ;cur_surResponse.BRMz.mean] ;
%     DELinst.BROop.mean = [instVec.BROop.mean ;cur_surResponse.BROop.mean] ;
%     DELinst.BRIp.mean  = [instVec.BRIp.mean ;cur_surResponse.BRIp.mean];
%     DELinst.TTMx.mean  = [instVec.TTMx.mean ;cur_surResponse.TTMx.mean] ;
%     DELinst.TTMy.mean  = [instVec.TTMy.mean ;cur_surResponse.TTMy.mean] ;
%     DELinst.TTMz.mean  = [instVec.TTMz.mean ;cur_surResponse.TTMz.mean] ;
%     DELinst.LSSMy.mean = [instVec.LSSMy.mean ;cur_surResponse.LSSMy.mean] ;
%     DELinst.LSSMz.mean = [instVec.LSSMz.mean ;cur_surResponse.LSSMz.mean] ;
%     DELinst.LSSTq.mean = [instVec.LSSTq.mean ;cur_surResponse.LSSTq.mean] ;
%
%     DELinst.TBMx.std  = [instVec.TBMx.std ;cur_surResponse.TBMx.std] ;
%     DELinst.TBMy.std  = [instVec.TBMy.std ;cur_surResponse.TBMy.std] ;
%     DELinst.TBMz.std  = [instVec.TBMz.std ;cur_surResponse.TBMz.std] ;
%     DELinst.BRMx.std  = [instVec.BRMx.std ;cur_surResponse.BRMx.std] ;
%     DELinst.BRMy.std  = [instVec.BRMy.std ;cur_surResponse.BRMy.std] ;
%     DELinst.BRMz.std  = [instVec.BRMz.std ;cur_surResponse.BRMz.std] ;
%     DELinst.BROop.std = [instVec.BROop.std ;cur_surResponse.BROop.std] ;
%     DELinst.BRIp.std  = [instVec.BRIp.std ;cur_surResponse.BRIp.std];
%     DELinst.TTMx.std  = [instVec.TTMx.std ;cur_surResponse.TTMx.std] ;
%     DELinst.TTMy.std  = [instVec.TTMy.std ;cur_surResponse.TTMy.std] ;
%     DELinst.TTMz.std  = [instVec.TTMz.std ;cur_surResponse.TTMz.std] ;
%     DELinst.LSSMy.std = [instVec.LSSMy.std ;cur_surResponse.LSSMy.std] ;
%     DELinst.LSSMz.std = [instVec.LSSMz.std ;cur_surResponse.LSSMz.std] ;
%     DELinst.LSSTq.std = [instVec.LSSTq.std ;cur_surResponse.LSSTq.std] ;
% end
% end

%
% function DAMinst = assignInstDAM(cur_surResponse,instVec,flagoff)
% if flagoff ==1
%     cur_surResponse.TBMx.mean = 0;
%     cur_surResponse.TBMy.mean = 0;
%     cur_surResponse.TBMz.mean = 0;
%     cur_surResponse.BRMx.mean = 0;
%     cur_surResponse.BRMy.mean = 0;
%     cur_surResponse.BRMz.mean = 0;
%     cur_surResponse.BROop.mean = 0;
%     cur_surResponse.BRIp.mean = 0;
%     cur_surResponse.TTMx.mean = 0;
%     cur_surResponse.TTMy.mean = 0;
%     cur_surResponse.TTMz.mean = 0;
%     cur_surResponse.TTMz.mean = 0;
%     cur_surResponse.LSSMy.mean = 0;
%     cur_surResponse.LSSMz.mean = 0;
%     cur_surResponse.LSSTq.mean = 0;
%
%     cur_surResponse.TBMx.std = 0;
%     cur_surResponse.TBMy.std = 0;
%     cur_surResponse.TBMz.std = 0;
%     cur_surResponse.BRMx.std = 0;
%     cur_surResponse.BRMy.std = 0;
%     cur_surResponse.BRMz.std = 0;
%     cur_surResponse.BROop.std = 0;
%     cur_surResponse.BRIp.std = 0;
%     cur_surResponse.TTMx.std = 0;
%     cur_surResponse.TTMy.std = 0;
%     cur_surResponse.TTMz.std = 0;
%     cur_surResponse.TTMz.std = 0;
%     cur_surResponse.LSSMy.std = 0;
%     cur_surResponse.LSSMz.std = 0;
%     cur_surResponse.LSSTq.std = 0;
% end
%
% if isempty(instVec)
%     DAMinst.TBMx.mean  = (cur_surResponse.TBMx.mean^4)*3600 ;
%     DAMinst.TBMy.mean  = (cur_surResponse.TBMy.mean^4)*3600 ;
%     DAMinst.TBMz.mean  = (cur_surResponse.TBMz.mean^4)*3600 ;
%     DAMinst.BRMx.mean  = (cur_surResponse.BRMx.mean^10)*3600 ;
%     DAMinst.BRMy.mean  = (cur_surResponse.BRMy.mean^10)*3600 ;
%     DAMinst.BRMz.mean  = (cur_surResponse.BRMz.mean^10)*3600 ;
%     DAMinst.BROop.mean = (cur_surResponse.BROop.mean^10)*3600 ;
%     DAMinst.BRIp.mean  = (cur_surResponse.BRIp.mean^10)*3600 ;
%     DAMinst.TTMx.mean  = (cur_surResponse.TTMx.mean^4)*3600 ;
%     DAMinst.TTMy.mean  = (cur_surResponse.TTMy.mean^4)*3600 ;
%     DAMinst.TTMz.mean  = (cur_surResponse.TTMz.mean^4)*3600 ;
%     DAMinst.LSSMy.mean = (cur_surResponse.LSSMy.mean^4)*3600 ;
%     DAMinst.LSSMz.mean = (cur_surResponse.LSSMz.mean^4)*3600 ;
%     DAMinst.LSSTq.mean = (cur_surResponse.LSSTq.mean^4)*3600 ;
%
%     DAMinst.TBMx.std  = (cur_surResponse.TBMx.std^4)*3600 ;
%     DAMinst.TBMy.std  = (cur_surResponse.TBMy.std^4)*3600 ;
%     DAMinst.TBMz.std  = (cur_surResponse.TBMz.std^4)*3600 ;
%     DAMinst.BRMx.std  = (cur_surResponse.BRMx.std^10)*3600 ;
%     DAMinst.BRMy.std  = (cur_surResponse.BRMy.std^10)*3600 ;
%     DAMinst.BRMz.std  = (cur_surResponse.BRMz.std^10)*3600 ;
%     DAMinst.BROop.std = (cur_surResponse.BROop.std^10)*3600 ;
%     DAMinst.BRIp.std  = (cur_surResponse.BRIp.std^10)*3600 ;
%     DAMinst.TTMx.std  = (cur_surResponse.TTMx.std^4)*3600 ;
%     DAMinst.TTMy.std  = (cur_surResponse.TTMy.std^4)*3600 ;
%     DAMinst.TTMz.std  = (cur_surResponse.TTMz.std^4)*3600 ;
%     DAMinst.LSSMy.std = (cur_surResponse.LSSMy.std^4)*3600 ;
%     DAMinst.LSSMz.std = (cur_surResponse.LSSMz.std^4)*3600 ;
%     DAMinst.LSSTq.std = (cur_surResponse.LSSTq.std^4)*3600 ;
% else
%     DAMinst.TBMx.mean  = [instVec.TBMx.mean ;(cur_surResponse.TBMx.mean^4)*3600] ;
%     DAMinst.TBMy.mean  = [instVec.TBMy.mean ;(cur_surResponse.TBMy.mean^4)*3600] ;
%     DAMinst.TBMz.mean  = [instVec.TBMz.mean ;(cur_surResponse.TBMx.mean^4)*3600];
%     DAMinst.BRMx.mean  = [instVec.BRMx.mean ;(cur_surResponse.BRMx.mean^10)*3600];
%     DAMinst.BRMy.mean  = [instVec.BRMy.mean ;(cur_surResponse.BRMy.mean^10)*3600];
%     DAMinst.BRMz.mean  = [instVec.BRMz.mean ;(cur_surResponse.BRMz.mean^10)*3600];
%     DAMinst.BROop.mean = [instVec.BROop.mean ;(cur_surResponse.BROop.mean^10)*3600];
%     DAMinst.BRIp.mean  = [instVec.BRIp.mean ;(cur_surResponse.BRIp.mean^10)*3600];
%     DAMinst.TTMx.mean  = [instVec.TTMx.mean ;(cur_surResponse.TTMx.mean^4)*3600] ;
%     DAMinst.TTMy.mean  = [instVec.TTMy.mean ;(cur_surResponse.TTMy.mean^4)*3600] ;
%     DAMinst.TTMz.mean  = [instVec.TTMz.mean ;(cur_surResponse.TTMz.mean^4)*3600] ;
%     DAMinst.LSSMy.mean = [instVec.LSSMy.mean ;(cur_surResponse.LSSMy.mean^4)*3600] ;
%     DAMinst.LSSMz.mean = [instVec.LSSMz.mean ;(cur_surResponse.LSSMz.mean^4)*3600] ;
%     DAMinst.LSSTq.mean = [instVec.LSSTq.mean ;(cur_surResponse.LSSTq.mean^4)*3600] ;
%
%     DAMinst.TBMx.std  = [instVec.TBMx.std ;(cur_surResponse.TBMx.std^4)*3600] ;
%     DAMinst.TBMy.std  = [instVec.TBMy.std ;(cur_surResponse.TBMy.std^4)*3600] ;
%     DAMinst.TBMz.std  = [instVec.TBMz.std ;(cur_surResponse.TBMx.std^4)*3600];
%     DAMinst.BRMx.std  = [instVec.BRMx.std ;(cur_surResponse.BRMx.std^10)*3600];
%     DAMinst.BRMy.std  = [instVec.BRMy.std ;(cur_surResponse.BRMy.std^10)*3600];
%     DAMinst.BRMz.std  = [instVec.BRMz.std ;(cur_surResponse.BRMz.std^10)*3600];
%     DAMinst.BROop.std = [instVec.BROop.std ;(cur_surResponse.BROop.std^10)*3600];
%     DAMinst.BRIp.std  = [instVec.BRIp.std ;(cur_surResponse.BRIp.std^10)*3600];
%     DAMinst.TTMx.std  = [instVec.TTMx.std ;(cur_surResponse.TTMx.std^4)*3600] ;
%     DAMinst.TTMy.std  = [instVec.TTMy.std ;(cur_surResponse.TTMy.std^4)*3600] ;
%     DAMinst.TTMz.std  = [instVec.TTMz.std ;(cur_surResponse.TTMz.std^4)*3600] ;
%     DAMinst.LSSMy.std = [instVec.LSSMy.std ;(cur_surResponse.LSSMy.std^4)*3600] ;
%     DAMinst.LSSMz.std = [instVec.LSSMz.std ;(cur_surResponse.LSSMz.std^4)*3600] ;
%     DAMinst.LSSTq.std = [instVec.LSSTq.std ;(cur_surResponse.LSSTq.std^4)*3600] ;
% end
% end
%
% function DAMcum = assignCumDEL(instVal,cumVal_prev)
%
% DAMcum.TBMx.mean  = [cumVal_prev.TBMx.mean; instVal.TBMx.mean(end)+ cumVal_prev.TBMx.mean(end)];
% DAMcum.TBMy.mean  = [cumVal_prev.TBMy.mean; instVal.TBMy.mean(end)+ cumVal_prev.TBMy.mean(end)];
% DAMcum.TBMz.mean  = [cumVal_prev.TBMz.mean; instVal.TBMz.mean(end)+ cumVal_prev.TBMz.mean(end)];
% DAMcum.BRMx.mean  = [cumVal_prev.BRMx.mean; instVal.BRMx.mean(end)+ cumVal_prev.BRMx.mean(end)];
% DAMcum.BRMy.mean  = [cumVal_prev.BRMy.mean; instVal.BRMy.mean(end)+ cumVal_prev.BRMy.mean(end)];
% DAMcum.BRMz.mean  = [cumVal_prev.BRMz.mean; instVal.BRMz.mean(end)+ cumVal_prev.BRMz.mean(end)];
% DAMcum.BROop.mean = [cumVal_prev.BROop.mean; instVal.BROop.mean(end)+ cumVal_prev.BROop.mean(end)];
% DAMcum.BRIp.mean  = [cumVal_prev.BRIp.mean; instVal.BRIp.mean(end)+ cumVal_prev.BRIp.mean(end)];
% DAMcum.TTMx.mean  = [cumVal_prev.TTMx.mean; instVal.TTMx.mean(end)+ cumVal_prev.TTMx.mean(end)];
% DAMcum.TTMy.mean  = [cumVal_prev.TTMy.mean; instVal.TTMy.mean(end)+ cumVal_prev.TTMy.mean(end)];
% DAMcum.TTMz.mean  = [cumVal_prev.TTMz.mean; instVal.TTMz.mean(end)+ cumVal_prev.TTMz.mean(end)];
% DAMcum.LSSMy.mean = [cumVal_prev.LSSMy.mean; instVal.LSSMy.mean(end)+ cumVal_prev.LSSMy.mean(end)];
% DAMcum.LSSMz.mean = [cumVal_prev.LSSMz.mean; instVal.LSSMz.mean(end)+ cumVal_prev.LSSMz.mean(end)];
% DAMcum.LSSTq.mean = [cumVal_prev.LSSTq.mean; instVal.LSSTq.mean(end)+ cumVal_prev.LSSTq.mean(end)];
% end
%
% function ControlActivity = getCntPercentage(CnttrCell)
% for i = 1:length(CnttrCell)
%     if strcmp(CnttrCell{i},'ShutDown')
%         countCont(i) = 0;
%     elseif strcmp(CnttrCell{i},'constTSR')
%         countCont(i) = 1;
%     elseif strcmp(CnttrCell{i},'constTSR_IPC')
%         countCont(i) = 2;
%     elseif strcmp(CnttrCell{i},'lin70')
%         countCont(i) = 3;
%     elseif strcmp(CnttrCell,'lin70_IPC')
%         countCont(i) = 4; %#ok<*AGROW>
%     end
% end
% ControlActivity.Shuts    = 100*sum(countCont==0)/length(CnttrCell);
% ControlActivity.constTSR = 100*sum(countCont==1)/length(CnttrCell);
% ControlActivity.constTSR_IPC = 100*sum(countCont==2)/length(CnttrCell);
% ControlActivity.li70     = 100*sum(countCont==3)/length(CnttrCell);
% ControlActivity.li70_IPC = 100*sum(countCont==4)/length(CnttrCell);
%
% end
