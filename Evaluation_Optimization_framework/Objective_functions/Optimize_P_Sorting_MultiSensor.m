function [Prat_out,optoptions,wohler] = Optimize_P_Sorting_MultiSensor (SortBins,BaseLine,AllBIns,surrogate,TSData,surrogate_case,binflag,FunSettings)

% Optimization function it should be in a standardize format and then copy
% paste to create different optimizers. The structure should be the same:
%
% - If bins find the most significant bins for DAM and REV
% - If all values just throw them all in
% - Assign actions (constraints to the output) to the xx most influential bins. Also according to being above or below rated
% - Objective function should come at the end of each function (as function in a function). This will make it easier to keep track
% - Feed to fmincon inputs, actions, options
% - when considering TS input instead the sorting is done only on below/above rated wind speed

for i =1:length(FunSettings.loadsens)
    if any(strcmp (FunSettings.loadsens{i},{'BROop';'BRIp';'BRMy';'BRMx';'BRMz'}))
        wohler(i) = 10;
    else
        wohler(i) = 4;
    end
    BaseDAM(i,1)  = BaseLine.cum.DAM.(FunSettings.loadsens{i});
end
BaseRev  = BaseLine.cum.Rev;

% Speed up
if surrogate_case==1
    [X1,Y1,Z1] = ndgrid(surrogate{1,1}.DataCnt.Dimensions.Dim1{2},surrogate{1, 1}.DataCnt.Dimensions.Dim2{2},surrogate{1,1}.DataCnt.Dimensions.Dim3{2});
end

%% Sort bins to revenue, dam and both and according to wsp
if binflag ==1

    IndOpt.Opt.All = SortBins.Rev_Base_sort(1:max(find(SortBins.Rev_Base_sort(:,1)>0)),2); %#ok<*MXFND> Find the non zero contributors. They are sorted so we need all indices until the last non zero
    % Find the indices that contribute the x% of each objective
    diffy.RevSigInd = find( cumsum(SortBins.Rev_Base_sort(:,1)/sum(SortBins.Rev_Base_sort(:,1))) >FunSettings.percToTotal, 1,"first" ); % from this indeces and up the contribution is higher than the percentage requested
    for i =1:length(FunSettings.loadsens)
        diffy.DamSigInd(i) = find( cumsum(SortBins.DAMBase.([FunSettings.loadsens{i} '_sort'])(:,1)/sum(SortBins.DAMBase.([FunSettings.loadsens{i} '_sort'])(:,1))) >FunSettings.percToTotal, 1,"first" );
    end
    % Separate to most influential in revenue and damage
    diffy.RevBins = SortBins.Rev_Base_sort(1:diffy.RevSigInd,2);
    diffy.DamBins =[];
    for i =1:length(FunSettings.loadsens)
        diffy.DamBinsInt = [diffy.DamBins;SortBins.DAMBase.([FunSettings.loadsens{i} '_sort'])(1:diffy.DamSigInd,2)];
    end
    diffy.DamBins = unique( diffy.DamBinsInt);
    RevOnly = setdiff(diffy.RevBins,diffy.DamBins);
    DamOnly = setdiff(diffy.DamBins,diffy.RevBins);

    % go through all points and assign actions and initial values according to wsp and objective
    for iBin = 1:length(IndOpt.Opt.All)
        curbin = IndOpt.Opt.All(iBin);
        curV = AllBIns.Vdisc(curbin);
        if curV>=11
            if ismember(curbin,RevOnly)
                % Assign bounds for x (limiting input variable space)
                lb(iBin) = 9.99;
                ub(iBin) = 13.01;
                % Assign initial values (rated values) to each bin
                x0(iBin) = 12;
            elseif ismember(curbin,DamOnly)
                lb(iBin) = 4.99;
                ub(iBin) = 10.01;
                x0(iBin) = 9;
            else
                lb(iBin) = 4.99;
                ub(iBin) = 13.01;
                x0(iBin) = 10;
            end
        else
            if ismember(curbin,RevOnly)
                lb(iBin) = 8.99;
                ub(iBin) = 10.01;
                x0(iBin) = 10;
            elseif ismember(curbin,DamOnly)
                lb(iBin) = 4.99;
                ub(iBin) = 10;
                x0(iBin) = 9;
            else
                lb(iBin) = 4.99;
                ub(iBin) = 10.01;
                x0(iBin) = 10;
            end
        end
    end
else
    IndOn = find(TSData.Price>0 & TSData.V>4 & TSData.V<24 & ~isnan(TSData.Price) )';
    IndOpt.Opt.All = IndOn;
    for iVal = 1:length(IndOpt.Opt.All)
        curV = TSData.V(IndOn(iVal));
        if curV>=11
            lb(iVal) = 4.99;
            ub(iVal) = 13;
            x0(iVal) = 10;
        else
            lb(iVal) = 4.99;
            ub(iVal) = 10.1;
            x0(iVal) = 10;
        end
    end
end

%% Apply Optimization

A = [];  % constrain the input according to actions for the region (A*x <= b)
b = [];  % constrain the actions according to actions for the region (A*x <= b)

Aeq = [];
beq = [];

options = optimoptions('fmincon', ...
    'DiffMinChange',0.01,...  % default 0
    'Algorithm','interior-point',...   %default  interior-point  / sqp /active-set
    'Diagnostics','on',...
    ...%     'FunValCheck','on',...
    'EnableFeasibilityMode',true,... % default false
    ...%     'MaxFunctionEvaluations',1e7, ...
    ...%     'FiniteDifferenceStepSize',0.01,... %default sqrt(eps) for forward finite differences, and eps^(1/3)
    'FiniteDifferenceType','forward',... % default  'central' 'forward'
    'Display','iter-detailed', ...  % 'iter' 'none'
    'FunValCheck','off', ...
    'MaxIterations',100,... %default 400; for the interior-point 600
    'StepTolerance',1e-14, ... % default 1e-6 for all except interior point 1e-10
    'MaxFunctionEvaluations',35000, ...% default 3000
    'UseParallel',true,...
    'ConstraintTolerance',1e-9... % default 1e-6
    ...%     'OptimalityTolerance',1e-12...  % default 1e-9
    );

optionsga = optimoptions('ga',....
    'Display','iter', ...  % 'iter'   'off' 'diagnose'
    'UseParallel',true,...
    'FunctionTolerance', 1.0e-03,... %default 1e-6
    'MaxGenerations',120, ...150*length(x0), ...  % default 100*length(inputs)
    'PopulationSize',1000, ...   % defaut 200
    'CrossoverFraction',0.8, ... % default 0.8 how many pass to the next gen
    'EliteCount',ceil(0.075*1000),... %default {ceil(0.05*PopulationSize)}
    'MigrationFraction', 0.3 ... % default 0.2
    ...%     'ConstraintTolerance',1e-1,... % default 1e-6
    ...%     'OptimalityTolerance',1e-12...  % default 1e-9
    );

optionsgam = optimoptions('gamultiobj',....
    'Display','iter', ...  % 'iter'   'off' 'diagnose'
    'UseParallel',true,...
    'FunctionTolerance', 1.0e-04,... %default 1e-6
    'MaxGenerations',1500, ...  % default 100*length(inputs)
    'PopulationSize',400, ...   % defaut 200
    'CrossoverFraction',0.8, ... % default 0.8 how many pass to the next gen
    'MigrationFraction', 0.2, ... % default 0.2
    'PlotFcn',{@gaplotpareto} ...
    );

if surrogate_case==1
    fun = @(x)func_obj(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.wDAM,FunSettings.wREV,IndOpt.Opt.All,FunSettings.loadsens,wohler,FunSettings.PenaltyRev,FunSettings.PenaltyDam,X1,Y1,Z1,binflag);
else
    fun = @(x)func_obj(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.wDAM,FunSettings.wREV,IndOpt.Opt.All,FunSettings.loadsens,wohler,FunSettings.PenaltyRev,FunSettings.PenaltyDam,[],[],[],binflag);
end

if surrogate_case==1
    funMulti = @(x)func_objMulti(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,IndOpt.Opt.All,FunSettings.loadsens,wohler,X1,Y1,Z1,binflag);
else
    funMulti = @(x)func_objMulti(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,IndOpt.Opt.All,FunSettings.loadsens,wohler,[],[],[],binflag);
end

%% Choose optimizer

if strcmp(FunSettings.Method, 'fmincon')
    %%%     FMINCON       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [x,~,~,output,~]  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
    optoptions.Options = options;
    toc

elseif strcmp(FunSettings.Method, 'GlobalSearch')
    gs = GlobalSearch( ...
        'NumTrialPoints', 1000, ...   default 1000
        'FunctionTolerance', 1e-3, ...  default 1e-6
        'BasinRadiusFactor' ,0.2, ... default 0.2
        'Display', 'iter', ...
        'NumTrialPoints', 10, ... default 1000
        'NumStageOnePoints', 10, ... default 200
        'XTolerance',  1e-4 ... default 1e-6
        );

    problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
    optoptions.Options = options;
    optoptions.gsOptions = gs;
    x = run(gs,problem);


elseif strcmp(FunSettings.Method, 'MultiStart')
    ms = MultiStart( ...
        'FunctionTolerance', 1e-3, ...  default 1e-6
        'Display', 'iter', ...
        'StartPointsToRun','all', ... default 'all'/ 'bounds'
        'UseParallel', false ...
        );
    problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
    optoptions.Options = options;
    [x,fval,exitflag,output,solutions] = run(ms,problem,5);

elseif strcmp(FunSettings.Method, 'GA')
    % %%% GENETIC ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [x,~,~,output,pops,scores]  = ga(fun,length(x0),A,b,Aeq,beq,lb,ub,[],optionsga); %#ok<*ASGLU>
    optoptions.Options = optionsga;
    toc

elseif strcmp(FunSettings.Method, 'GAmulti')
    %%% MULTI OBJECTIVE GA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [xout,fval,~,output] = gamultiobj(funMulti,length(x0),A,b,Aeq,beq,lb,ub,[],optionsgam);
    optoptions.Options = optionsgam;
    optoptions.Pareto = fval;
    figure,plot(-fval(:,1),fval(:,2),'o'),ylabel('Load Delta to baseline %'),xlabel('Revenue Delta to baseline %'),grid on
    allval = find(fval(:,1)<0 & fval(:,2)<0);
    if ~isempty(allval)
        x = xout(find(fval(:,1)==min(fval(allval,1))),:)'; %#ok<*FNDSB>
    else
        x = xout(find(fval(:,1)==min(fval(:,1))),:)';
    end
    toc

end

%% unpacking the bins to TS and assign outputs
x = round(x*10)/10;
Prat_out = zeros(length(TSData.V),1);
if binflag==1
    for iT = 1:length(TSData.V)
        if TSData.Price(iT)<=0 || TSData.V(iT)<=4 || TSData.V(iT)>=24 || isnan(TSData.Price(iT))
            Prat_out(iT,1) =0;
        else
            row =  discretize(TSData.V(iT),AllBIns.bin.V.edge) ;
            if isnan(row)
                row = length(AllBIns.bin.V.center);
            end
            col =  discretize(TSData.Price(iT),AllBIns.bin.Price.edge) ;
            Prat_out(iT,1) = x(IndOpt.Opt.All==sub2ind([length(AllBIns.bin.V.center) length(AllBIns.bin.Price.center)],row,col));
        end
    end
else
    Prat_out(IndOn) = x;
end

if exist("output","var")
    optoptions.Output = output;
end
optoptions.Values = x;

end

%% Objective function

function obj = func_obj(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,wDAM,wRev,BinInd,loadsens,wohler,PenaltyRev,PenaltyDam,X1,Y1,Z1,binflag)

% 1. Loop over time and assign a power to each hour block
% 2. Probe once the surrogate to take all the responses
% 3. Calculate cumulative Damage and Revenue
% 4. Calculate Deltas of REV and DAM compared to the baseline
% 5. Combine the two deltas with weight in one
% 6. Check for surpassing thresholds and add penalty if they do

IndShut = find(TSData.Price<=0 | TSData.V<4 | TSData.V>24 | isnan(TSData.Price) )';
IndOn = setdiff(1:length(TSData.V),IndShut);
P = zeros(length(TSData.V),1);
if binflag==1
    rows = discretize(TSData.V(IndOn),AllBIns.bin.V.edge) ;
    rows(isnan(rows)) = length(AllBIns.bin.V.center);
    cols =  discretize(TSData.Price(IndOn),AllBIns.bin.Price.edge) ;
    OpInd_bins =sub2ind([length(AllBIns.bin.V.center) length(AllBIns.bin.Price.center)],rows,cols); % Which value of the TS is in which bin for all timesteps by index
    for i=1:length(IndOn)
        P(IndOn(i),1) = x(BinInd==OpInd_bins(i)); %#ok<*AGROW,*NASGU> % assign the poewr level according to the bins
    end
else
    P(IndOn) = x;
end
if surrogate_case == 1
    SurOut.Power.mean  = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.Power.mean,TSData.V(IndOn),TSData.TI(IndOn),P(IndOn),'spline'));
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean  = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.(loadsens{i}).mean,TSData.V(IndOn),TSData.TI(IndOn),P(IndOn),'spline'));
    end
elseif surrogate_case == 2
    SurOut.Power.mean = predict(surrogate{1,1}.Output.Power,[TSData.V(IndOn),TSData.TI(IndOn),P(IndOn)]);
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean  = predict(surrogate{1,1}.Output.(loadsens{i}),[TSData.V(IndOn),TSData.TI(IndOn),P(IndOn)]);
    end
end

% Calculate new cummulative values
RevCum = sum(TSData.Price(IndOn).*SurOut.Power.mean/1000);
for i =1:length(loadsens)
    DamCum(i,1) = sum(3600*(SurOut.(loadsens{i}).mean).^wohler(i));
end
%Calculate the relative differencs
REV  = 100*(RevCum/BaseRev-1);
if REV>1  % used for load optimization to restrict revenue I used 1 for revenue
    REV=REV;
end
DAMint = 100*(DamCum./BaseDAM-1);
dlim = -10;
IndHigh = find(DAMint<dlim);
DAMInt1 = DAMint;
if ~isempty(IndHigh)
    DAMint(IndHigh) = dlim; % to avoid load reductions dominating the single objective I used 10 for revenue and 20 for load
end
DAM  = mean(DAMint);   % here different types of weighting can be introduced

% Combining the objectives to one with weights
obj = -wRev*REV+wDAM*DAM;

if any(DAMInt1>PenaltyDam)
        IndHigh2 = find(DAMInt1>PenaltyDam);
        obj = obj + 10*length(IndHigh2);
end
if DAM > PenaltyDam
    %     obj = abs(obj*DAM); % use for revenue
    obj = obj+10; % use for load
    
elseif REV<-2*PenaltyRev
    obj = 10*abs(obj*(REV+2)); % use for revenue
% obj = obj+10*abs(1+abs(REV/(PenaltyRev+0.1))); % use for load
elseif REV<PenaltyRev
    obj = abs(obj*(REV+2)); % use for revenue
%     obj = obj+10*abs(1+abs(REV/(PenaltyRev+0.1)));  % use for load
end


end

function obj = func_objMulti(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,BinInd,loadsens,wohler,X1,Y1,Z1,binflag)

% 1. Loop over time and assign a power to each hour block
% 2. Probe once the surrogate to take all the responses
% 3. Calculate cumulative Damage and Revenue
% 4. Calculate Deltas of REV and DAM compared to the baseline
% 5. Combine the two deltas with weight in one
% 6. Check for surpassing thresholds and add penalty if they do

IndShut = find(TSData.Price<=0 | TSData.V<4 | TSData.V>24 | isnan(TSData.Price) )';
IndOn = setdiff(1:length(TSData.V),IndShut);
P = zeros(length(TSData.V),1);
if binflag==1
    rows = discretize(TSData.V(IndOn),AllBIns.bin.V.edge) ;
    rows(isnan(rows)) = length(AllBIns.bin.V.center);
    cols =  discretize(TSData.Price(IndOn),AllBIns.bin.Price.edge) ;
    OpInd_bins =sub2ind([length(AllBIns.bin.V.center) length(AllBIns.bin.Price.center)],rows,cols); % Which value of the TS is in which bin for all timesteps by index
    for i=1:length(IndOn)
        P(IndOn(i),1) = x(BinInd==OpInd_bins(i)); %#ok<*AGROW,*NASGU> % assign the poewr level according to the bins
    end
else
    P(IndOn) = x;
end
if surrogate_case == 1
    SurOut.Power.mean  = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.Power.mean,TSData.V(IndOn),TSData.TI(IndOn),P(IndOn),'spline'));
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean  = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.(loadsens{i}).mean,TSData.V(IndOn),TSData.TI(IndOn),P(IndOn),'spline'));
    end
elseif surrogate_case == 2
    SurOut.Power.mean = predict(surrogate{1,1}.Output.Power,[TSData.V(IndOn),TSData.TI(IndOn),P(IndOn)]);
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean  = predict(surrogate{1,1}.Output.(loadsens{i}),[TSData.V(IndOn),TSData.TI(IndOn),P(IndOn)]);
    end
end

% Calculate new cummulative values
RevCum = sum(TSData.Price(IndOn).*SurOut.Power.mean/1000);
for i =1:length(loadsens)
    DamCum(i,1) = sum(3600*(SurOut.(loadsens{i}).mean).^wohler(i));
end
%Calculate the relative differencs
REV  = 100*(RevCum/BaseRev-1);
DAMint = 100*(DamCum./BaseDAM-1);
obj(1)= -REV;
for i=1:length(DAMint)
    obj(end+1)= DAMint(i);
end

end

