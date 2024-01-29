function [Prat_out,optoptions,wohler] = Optimize_P_Shut_IBC_Sorting_MultiSensor (SortBins,BaseLine,AllBIns,surrogate,TSData,surrogate_case,binflag,FunSettings)

% Optimization function it should be in a standardize format and then copy
% paste to create different optimizers. The structure should be the same:
%
% - If bins find the most significant bins for DAM and REV
% - If all values just throw them all in
% - Assign actions (constraints to the output) to the xx most influential bins
% - Objective function should come at the end of each function (as function in a function). This will make it easier to keep track
% - Feed to fmincon inputs, actions, options
% - when considering TS input instead the sorting is done only on below/above rated wind speed

NoOpVals = 0;
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

%% Prepare inputs for objective function to speed up

% Identify the indices of the operating points
Ind.Shut = find(TSData.Price<=0 | TSData.V<=4 | TSData.V>=24 | isnan(TSData.Price) )';
Ind.On = setdiff(1:length(TSData.V),Ind.Shut);

% sort bins or values according to their contribution to the obkectives and assign relevant actions
if binflag ==1
    IndOpt.Opt.All = SortBins.Rev_Base_sort(1:max(find(SortBins.Rev_Base_sort(:,1)>0)),2); %#ok<*MXFND> Find the non zero contributors

    % Find the indices that contribute the x% of each objective
    diffy.RevSigInd = find( cumsum(SortBins.Rev_Base_sort(:,1)/sum(SortBins.Rev_Base_sort(:,1))) >FunSettings.percToTotal, 1,"first" );
    for i =1:length(FunSettings.loadsens)
        diffy.DamSigInd = find( cumsum(SortBins.DAMBase.([FunSettings.loadsens{i} '_sort'])(:,1)/sum(SortBins.DAMBase.([FunSettings.loadsens{i} '_sort'])(:,1))) >FunSettings.percToTotal, 1,"first" );
    end
    % Separate to most influential in revenue and damage
    diffy.RevBins = SortBins.Rev_Base_sort(1:diffy.RevSigInd,2);
    for i =1:length(FunSettings.loadsens)
        diffy.DamBins = SortBins.DAMBase.([FunSettings.loadsens{i} '_sort'])(1:diffy.DamSigInd,2);
    end
    RevOnly = setdiff(diffy.RevBins,diffy.DamBins);
    DamOnly = setdiff(diffy.DamBins,diffy.RevBins);

    % go through all points and assign actions and initial values according to wsp and objective
    for iBin = 1:length(IndOpt.Opt.All)
        curbin = IndOpt.Opt.All(iBin);
        curV = AllBIns.Vdisc(curbin);
        if curV>=11
            if ismember(curbin,RevOnly)
                % Assign bounds for x (limiting input variable space)
                lb(iBin) = 9.0;
                ub(iBin) = 13.01;
                % Assign initial values (rated values) to each bin
                x0(iBin) = 11;
            elseif ismember(curbin,DamOnly)
                lb(iBin) = 4.99;
                ub(iBin) = 11.51;
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
    if FunSettings.constval.Shut==1
    else
        x0(end+1) = FunSettings.constval.ShutVal;  % input for price thrweshold to shut down
        lb(end+1) = 0; % eur revenue threshold  (also a TI theshold here)
        ub(end+1) = 40; % eur revenue threshold  (also a TI theshold here)
    end
    if FunSettings.constval.IBC==1
    else
        x0(end+1) = FunSettings.constval.IBCVal;  % input for IBC to start activating (it will not work with this sensor...)
        lb(end+1) = 13; %  m/s wsp to activate IBC (I could also add a TI threshold here)
        ub(end+1) = 24;  %  m/s wsp to activate IBC (I could also add a TI threshold here)
    end
else
    IndOn = find(TSData.Price>0 & TSData.V>4 & TSData.V<24 & ~isnan(TSData.Price) )';
    IndOpt.Opt.All = IndOn;
    if ~isempty(IndOn)
        for iVal = 1:length(IndOpt.Opt.All)
            curV = TSData.V(IndOn(iVal));
            if curV>=11
                lb(iVal) = 4.99;
                ub(iVal) = 13.01;
                x0(iVal) = 10;
            else
                lb(iVal) = 5;
                ub(iVal) = 10;
                x0(iVal) = 10;
            end
        end
        if FunSettings.constval.Shut==1
        else
            x0(end+1) = FunSettings.constval.ShutVal;  % input for price thrweshold to shut down
            lb(end+1) = 0; % eur revenue threshold  (also a TI theshold here)
            ub(end+1) = 40; % eur revenue threshold  (also a TI theshold here)
        end
        if FunSettings.constval.IBC==1
        else
            x0(end+1) = FunSettings.constval.IBCVal;  % input for IBC to start activating (it will not work with this sensor...)
            lb(end+1) = 13; %  m/s wsp to activate IBC (I could also add a TI threshold here)
            ub(end+1) = 24;  %  m/s wsp to activate IBC (I could also add a TI threshold here)
        end
    else
        NoOpVals = 1;
    end
end

%% Set  Optimization settings

A = [];  % constrain the input according to actions for the region (A*x <= b)
b = [];  % constrain the actions according to actions for the region (A*x <= b)

Aeq = []; % linear equality constrains. Not applicable for me
beq = []; % linear equality constrains. Not applicable for me

options = optimoptions('fmincon', ...
    'DiffMinChange',0.1,...  % default 0
    'Algorithm','interior-point',...   %default  interior-point  / sqp /active-set
    'Diagnostics','on',...
    ...%     'FunValCheck','on',...
    'EnableFeasibilityMode',true,... % default false
    ...%     'MaxFunctionEvaluations',1e7, ...
    ...%     'FiniteDifferenceStepSize',0.01,... %default sqrt(eps) for forward finite differences, and eps^(1/3)
    'FiniteDifferenceType','forward',... % default  'central'
    'Display','iter-detailed', ...  % 'iter' 'none'
    'FunValCheck','off', ...
    'MaxIterations',600,... %default 400; for the interior-point 600
    'StepTolerance',1e-11, ... % default 1e-6 for all except interior point 1e-10
    'MaxFunctionEvaluations',250000, ...% default 3000
    'UseParallel',true...
    ...%     'ConstraintTolerance',1e-1,... % default 1e-6
    ...%     'OptimalityTolerance',1e-12...  % default 1e-9
    );

optionsga = optimoptions('ga',....
    'Display','iter', ...  % 'iter'   'off' 'diagnose'
    'UseParallel',true,...
    'FunctionTolerance', 1.0e-03,... %default 1e-6
    'MaxGenerations',150, ...  % default 100*length(inputs)
    'PopulationSize',1500, ...   % defaut 200  used 150
    'CrossoverFraction',0.95, ... % default 0.8 how many pass to the next gen
    'EliteCount',ceil(0.075*1500) ... %default {ceil(0.05*PopulationSize)}
    ... %     'MigrationFraction', 0.2 ... % default 0.2  used 0.3
    ...%     'ConstraintTolerance',1e-1,... % default 1e-6
    ...%     'OptimalityTolerance',1e-12...  % default 1e-9
    );

optionsPS = optimoptions("patternsearch", ...     % https://www.mathworks.com/help/gads/patternsearch.html#buxdit7-2
    'Display','iter', ...  % 'iter'   'off' 'diagnose'
    'FunctionTolerance', 1.0e-03,... %default 1e-6
    'MaxIterations',600,... %default 100*nvars;
    'UseParallel',true,...
    'MeshTolerance', 1e-4, ... % defailt 1e-6
    'PollMethod','GPSPositiveBasis2N', ... % default 'GPSPositiveBasis2N','GPSPositiveBasisNp1','GSSPositiveBasisNp1', 'MADSPositiveBasis2N', 'MADSPositiveBasisNp1'
    'StepTolerance', 1e-4,... %default 1e-6
    'MeshContractionFactor', 0.5,... %default 0.5
    'MeshExpansionFactor',1.2, ... % default 2
    'InitialMeshSize',0.2,... % default 1
    'MaxMeshSize',Inf,... % default Inf
    'ScaleMesh',true ... % default true
    );

optionsParticle = optimoptions("particleswarm", ...     % https://www.mathworks.com/help/gads/patternsearch.html#buxdit7-2
    'Display','iter', ...  % 'iter'   'off' 'diagnose'
    'FunctionTolerance',5e-04,... %default 1e-6
    'MaxIterations',250,... %default 200*nvars;
    'UseParallel',true,... %
    'MaxStallIterations',12,... % 20
    'SelfAdjustmentWeight', 1.49, ... % 1.49
    'SocialAdjustmentWeight',1.49,... % 1.49
    'SwarmSize',600, ... % 10*nvars
    'ObjectiveLimit',-35,...% -Inf
    'HybridFcn', []... % []  {@fmincon,options}
    );


optionsgam = optimoptions('gamultiobj',....
    'Display','iter', ...  % 'iter'   'off' 'diagnose'
    'UseParallel',true,...
    'FunctionTolerance', 1.0e-04,... %default 1e-6
    'MaxGenerations',1500, ...  % default 100*length(inputs)
    'PopulationSize',400, ...   % defaut 200
    'CrossoverFraction',0.75, ... % default 0.8 how many pass to the next gen
    'MigrationFraction', 0.2, ... % default 0.2
    'PlotFcn',{@gaplotpareto} ...
    );

if surrogate_case==1
    fun = @(x)func_obj2(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.wDAM,FunSettings.wREV,IndOpt.Opt.All,FunSettings.loadsens,wohler,FunSettings.PenaltyRev,FunSettings.PenaltyDam,X1,Y1,Z1,binflag,Ind,FunSettings);
else
    fun = @(x)func_obj2(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.wDAM,FunSettings.wREV,IndOpt.Opt.All,FunSettings.loadsens,wohler,FunSettings.PenaltyRev,FunSettings.PenaltyDam,[],[],[],binflag,Ind,FunSettings);
end

if surrogate_case==1
    funMulti = @(x)func_objMulti(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,IndOpt.Opt.All,FunSettings.loadsens,wohler,X1,Y1,Z1,binflag,Ind,FunSettings);
else
    funMulti = @(x)func_objMulti(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,IndOpt.Opt.All,FunSettings.loadsens,wohler,[],[],[],binflag,Ind,FunSettings);
end

%% Choose optimizer
if NoOpVals ==1
    x = 0*ones(length(TSData.V),1);
    if FunSettings.constval.Shut==1
        x(end+1) = 0;
    end
    if FunSettings.constval.IBC==1
        x(end+1) = 0;
    end
    fval =-2;
    optoptions.fval =fval;
elseif strcmp(FunSettings.Method, 'fmincon')
    %%%     FMINCON       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [x,fval,~,output,~]  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
    optoptions.Options = options;
    optoptions.fval =fval;
    toc

    %%%   GLOBAL SEARCH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(FunSettings.Method, 'GlobalSearch')
    gs = GlobalSearch( ...
        'NumTrialPoints', 1500, ...   default 1000
        'FunctionTolerance', 1e-3, ...  default 1e-6
        'BasinRadiusFactor' ,0.2, ... default 0.2
        'Display', 'iter', ...  'none' 'iter'
        'NumStageOnePoints', 400, ... default 200
        'XTolerance', 1e-4, ... default 1e-6
        'MaxTime',1500 ...
        );
    problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
    optoptions.Options = options;
    [x,fval] = run(gs,problem);
    optoptions.fval =fval;

elseif strcmp(FunSettings.Method, 'MultiStart')
    %%% MULTISTART    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ms = MultiStart( ...
        'FunctionTolerance', 1e-3, ...  default 1e-6
        'Display', 'iter', ...
        'StartPointsToRun','all', ... default 'all'/ 'bounds'
        'UseParallel', true ...
        );
    problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
    x = run(ms,problem,100);

elseif strcmp(FunSettings.Method, 'GA')
    %%% GENETIC ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [x,fval,~,output,pops,scores]  = ga(fun,length(x0),A,b,Aeq,beq,lb,ub,[],optionsga); %#ok<*ASGLU>
    optoptions.Options = optionsga;
    optoptions.fval =fval;
    toc

elseif strcmp(FunSettings.Method, 'PatternSearch')
    %%% PATTERN SEARCH    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [x,fval,~,output]  = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,[],optionsPS); %#ok<*ASGLU>
    optoptions.Options = optionsPS;
    optoptions.fval =fval;
    toc

elseif strcmp(FunSettings.Method, 'ParticleSwarm')
    %%% PARTICLE SWARM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [x,fval,~,output]  = particleswarm(fun,length(x0),lb,ub,optionsParticle); %#ok<*ASGLU>
    optoptions.Options = optionsParticle;
    optoptions.fval =fval;
    toc

elseif strcmp(FunSettings.Method, 'GAmulti')
    %%% MULTI OBJECTIVE GA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [xout,fval,~,output] = gamultiobj(funMulti,length(x0),A,b,Aeq,beq,lb,ub,[],optionsgam);
    optoptions.Options = optionsgam;
    optoptions.Pareto = fval;
    if exist(AllbIns,"var")
        optoptions.Bins = AllBIns.bin;
    end
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
    Prat_out = zeros(length(TSData.V),1);
    Prat_out(Ind.On) = x(1:length(Ind.On));
end

% Assign thresholds
if FunSettings.constval.Shut==0 && FunSettings.constval.IBC==0
    Thresholds.IBC.WSP = x(end);
    Thresholds.Shut.Rev = x(end-1);
elseif FunSettings.constval.IBC==1 && FunSettings.constval.Shut==1
    Thresholds.IBC.WSP = FunSettings.constval.IBCVal;
    Thresholds.Shut.Rev = FunSettings.constval.ShutVal;
elseif FunSettings.constval.IBC==0 && FunSettings.constval.Shut==1
    Thresholds.IBC.WSP = x(end);
    Thresholds.Shut.Rev = FunSettings.constval.ShutVal;
elseif FunSettings.constval.IBC==1 && FunSettings.constval.Shut==0
    Thresholds.IBC.WSP = FunSettings.constval.IBCVal;
    Thresholds.Shut.Rev = x(end);
end


if exist("output","var")
    optoptions.Output = output;
end
optoptions.Values = x;
optoptions.Thresholds = Thresholds;


end

%% Objective function
function obj = func_obj2(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,wDAM,wRev,BinInd,loadsens,wohler,PenaltyRev,PenaltyDam,X1,Y1,Z1,binflag,Ind,FunSettings)
% 1. Check IBC on/off according to thresholds for the whole TS
% 2. Assign the power level from the optimizer for thte whole TS unpacking the bins
% 3. Probe twice the surrogate to take all the responses of the TS for the noIBC and IBC time steps based on 1
% 4. Go throuth the time series and find the shut down points according to the shut thresholds
% 5.Change their response REV and DAM to 0 for these points
% 6. Calculate cumulative Damage and Revenue
% 7. Calculate Deltas of REV and DAM compared to the baseline
% 8. Combine the two deltas with weight in one
% 9. Check for surpassing thresholds and add penalty if they do

% separate the values to P and thresholds
if length(x)  == length(BinInd)
    IBC_thres = FunSettings.constval.IBCVal;
    Shut_thres = FunSettings.constval.ShutVal;
elseif length(x)  == length(BinInd)+2
    IBC_thres  = x(end);
    Shut_thres = x(end-1);
    x = x(1:end-2);
elseif FunSettings.constval.Shut==1
    Shut_thres = FunSettings.constval.ShutVal;
    IBC_thres = x(end);
    x = x(1:end-1);
elseif FunSettings.constval.IBC==1
    IBC_thres  = FunSettings.constval.IBCVal;
    Shut_thres = x(end);
    x = x(1:end-1);
end

% Assign IBC according to the new thresholds
Ind.newIBC = find(TSData.Price>0 & TSData.V>=IBC_thres & TSData.V<=24 & ~isnan(TSData.Price));
Ind.noIBC = setdiff(Ind.On,Ind.newIBC);
Ind.IBC   = setdiff(Ind.On,Ind.noIBC);

P = zeros(length(TSData.V),1); % all power rating of TS
if binflag==1
    % assign optimizer power levels
    rows = discretize(TSData.V(Ind.On),AllBIns.bin.V.edge) ;
    rows(isnan(rows)) = length(AllBIns.bin.V.center);
    cols =  discretize(TSData.Price(Ind.On),AllBIns.bin.Price.edge) ;
    OpInd_bins =sub2ind([length(AllBIns.bin.V.center) length(AllBIns.bin.Price.center)],rows,cols);% Which value of the time series is in which bin for all operational timesteps
    clear rows cols
    for i=1:length(Ind.On)
        P(Ind.On(i),1) = x(BinInd==OpInd_bins(i)); %#ok<*AGROW,*NASGU>
    end
else
    P(Ind.On) = x;
end
IndLow = find(TSData.V<=5.5 & P<10  );
if ~isempty(IndLow)
    P(IndLow) = 7.5;
end
% Get the surrpgate response for each time step for both cases
SurOut.Power.mean =  zeros(length(TSData.V),1);
for i =1:length(loadsens)
    SurOut.(loadsens{i}).mean = zeros(length(TSData.V),1);
end
if surrogate_case == 1
    % no IBC
    SurOut.Power.mean(Ind.noIBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.Power.mean,TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC),'spline'));
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean(Ind.noIBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.(loadsens{i}).mean,TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC),'spline'));
        % IBC
        SurOut.(loadsens{i}).mean(Ind.IBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,2}.DataCnt.(loadsens{i}).mean,TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC),'spline'));
    end
    SurOut.Power.mean(Ind.IBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,2}.DataCnt.Power.mean,TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC),'spline'));

elseif surrogate_case == 2
    % no IBC
    SurOut.Power.mean(Ind.noIBC) = predict(surrogate{1,1}.Output.Power,[TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC)]);
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean(Ind.noIBC) = predict(surrogate{1,1}.Output.(loadsens{i}),[TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC)]);
        % IBC
        SurOut.(loadsens{i}).mean(Ind.IBC) = predict(surrogate{1,2}.Output.(loadsens{i}),[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
    end
    SurOut.Power.mean(Ind.IBC) = predict(surrogate{1,2}.Output.Power,[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
end
if any(SurOut.Power.mean <0)
    SurOut.Power.mean(SurOut.Power.mean <0) = 10;
end
% Find shut down bins according to the new threshold and remove them from the list to be evaluated
RevProj = (SurOut.Power.mean/1000).*TSData.Price;
Ind.ShutFinal = find(RevProj<Shut_thres);
SurOut.Power.mean(Ind.ShutFinal) = 0;
for i =1:length(loadsens)
    SurOut.(loadsens{i}).mean(Ind.ShutFinal) = 0;
end

% Calculate new cummulative values
RevCum = sum(TSData.Price.*SurOut.Power.mean/1000);
for i =1:length(loadsens)
    DamCum(i,1) = sum(3600*(SurOut.(loadsens{i}).mean).^wohler(i));
end
%Calculate the relative differencs
REV  = 100*(RevCum/BaseRev-1);
revlim = FunSettings.revlim;
if REV>revlim  % used for load optimization to restrict revenue I used 1 for revenue
    REV=revlim;
end
DAMint = 100*(DamCum./BaseDAM-1);
dlim = FunSettings.dlim;
IndHigh = find(DAMint<dlim);
DAMInt1 = DAMint;
if ~isempty(IndHigh)
    DAMint(IndHigh) = dlim; % to avoid load reductions dominating the single objective I used 10 for revenue and 20 for load
end
DAM  = mean(DAMint);

% Combining the objectives to one with weights
obj = -wRev*REV+wDAM*DAM;

% Penalties to force the optimizer to some direction
if any(DAMInt1>PenaltyDam)
    IndHigh2 = find(DAMInt1>PenaltyDam);
    obj = obj + 10*length(IndHigh2);
end
if DAM > PenaltyDam
    %     obj = abs(obj*DAM); % use for revenue
    obj = obj+10; % use for load

end

if REV<PenaltyRev
    if REV<0
        obj = obj +100*(2+abs((abs(REV)+1)/(PenaltyRev+0.1)));%10*abs(obj*(REV+5)); % use for revenue
    else
        %         if REV>=1
        %         obj = obj +10*(2+abs((REV+1)/(PenaltyRev+0.1)));%10*abs(obj*(REV+5)); % use for revenue
        %         else
        obj = obj +10*(2+abs((10*REV+1)/(PenaltyRev+0.1)));%10*abs(obj*(REV+5)); % use for revenue
        %         end
    end
    % obj = obj+10*abs(1+abs(REV/(PenaltyRev+0.1))); % use for load
    % elseif REV<PenaltyRev
    %     obj = abs(obj*(REV+5)); % use for revenue
    %     obj = obj+10*abs(1+abs(REV/(PenaltyRev+0.1)));  % use for load
end

end


function  obj = func_objMulti(x,AllBIns,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,BinInd,loadsens,wohler,X1,Y1,Z1,binflag,Ind,FunSettings)
% 1. Check IBC on/off according to thresholds for the whole TS
% 2. Assign the power level from the optimizer for thte whole TS unpacking the bins
% 3. Probe twice the surrogate to take all the responses of the TS for the noIBC and IBC time steps based on 1
% 4. Go throuth the time series and find the shut down points according to the shut thresholds
% 5.Change their response REV and DAM to 0 for these points
% 6. Calculate cumulative Damage and Revenue
% 7. Calculate Deltas of REV and DAM compared to the baseline
% 8. Combine the two deltas with weight in one
% 9. Check for surpassing thresholds and add penalty if they do

% separate the values to P and thresholds
if length(x)  == length(BinInd)
    IBC_thres = FunSettings.constval.IBCVal;
    Shut_thres = FunSettings.constval.ShutVal;
elseif length(x)  == length(BinInd)+2
    IBC_thres  = x(end);
    Shut_thres = x(end-1);
    x = x(1:end-2);
elseif FunSettings.constval.Shut==1
    Shut_thres = FunSettings.constval.ShutVal;
    IBC_thres = x(end);
    x = x(1:end-1);
elseif FunSettings.constval.IBC==1
    IBC_thres  = FunSettings.constval.IBCVal;
    Shut_thres = x(end);
    x = x(1:end-1);
end

% Assign IBC according to the new thresholds
Ind.newIBC = find(TSData.Price>0 & TSData.V>=IBC_thres & TSData.V<=24 & ~isnan(TSData.Price));
Ind.noIBC = setdiff(Ind.On,Ind.newIBC);
Ind.IBC   = setdiff(Ind.On,Ind.noIBC);

P = zeros(length(TSData.V),1); % all power rating of TS
if binflag==1
    % assign optimizer power levels
    rows = discretize(TSData.V(Ind.On),AllBIns.bin.V.edge) ;
    rows(isnan(rows)) = length(AllBIns.bin.V.center);
    cols =  discretize(TSData.Price(Ind.On),AllBIns.bin.Price.edge) ;
    OpInd_bins =sub2ind([length(AllBIns.bin.V.center) length(AllBIns.bin.Price.center)],rows,cols);% Which value of the time series is in which bin for all operational timesteps
    clear rows cols
    for i=1:length(Ind.On)
        P(Ind.On(i),1) = x(BinInd==OpInd_bins(i)); %#ok<*AGROW,*NASGU>
    end
else
    P(Ind.On) = x;
end
IndLow = find(TSData.V<=5.5 & P<10  );
if ~isempty(IndLow)
    P(IndLow) = 7.5;
end
% Get the surrpgate response for each time step for both cases
SurOut.Power.mean =  zeros(length(TSData.V),1);
for i =1:length(loadsens)
    SurOut.(loadsens{i}).mean = zeros(length(TSData.V),1);
end
if surrogate_case == 1
    % no IBC
    SurOut.Power.mean(Ind.noIBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.Power.mean,TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC),'spline'));
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean(Ind.noIBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.(loadsens{i}).mean,TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC),'spline'));
        % IBC
        SurOut.(loadsens{i}).mean(Ind.IBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,2}.DataCnt.(loadsens{i}).mean,TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC),'spline'));
    end
    SurOut.Power.mean(Ind.IBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,2}.DataCnt.Power.mean,TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC),'spline'));

elseif surrogate_case == 2
    % no IBC
    SurOut.Power.mean(Ind.noIBC) = predict(surrogate{1,1}.Output.Power,[TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC)]);
    for i =1:length(loadsens)
        SurOut.(loadsens{i}).mean(Ind.noIBC) = predict(surrogate{1,1}.Output.(loadsens{i}),[TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC)]);
        % IBC
        SurOut.(loadsens{i}).mean(Ind.IBC) = predict(surrogate{1,2}.Output.(loadsens{i}),[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
    end
    SurOut.Power.mean(Ind.IBC) = predict(surrogate{1,2}.Output.Power,[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
end
if any(SurOut.Power.mean <0)
    SurOut.Power.mean(SurOut.Power.mean <0) = 10;
end
% Find shut down bins according to the new threshold and remove them from the list to be evaluated
RevProj = (SurOut.Power.mean/1000).*TSData.Price;
Ind.ShutFinal = find(RevProj<Shut_thres);
SurOut.Power.mean(Ind.ShutFinal) = 0;
for i =1:length(loadsens)
    SurOut.(loadsens{i}).mean(Ind.ShutFinal) = 0;
end

% Calculate new cummulative values
RevCum = sum(TSData.Price.*SurOut.Power.mean/1000);
for i =1:length(loadsens)
    DamCum(i,1) = sum(3600*(SurOut.(loadsens{i}).mean).^wohler(i));
end

%Calculate the relative differencs
REV  = 100*(RevCum/BaseRev-1);
revlim = FunSettings.revlim;
if REV>revlim  % used for load optimization to restrict revenue I used 1 for revenue
    REV=revlim;
end
DAMint = 100*(DamCum./BaseDAM-1);
dlim = FunSettings.dlim;
IndHigh = find(DAMint<dlim);
DAMInt1 = DAMint;
if ~isempty(IndHigh)
    DAMint(IndHigh) = dlim; % to avoid load reductions dominating the single objective I used 10 for revenue and 20 for load
end
DAM  = mean(DAMint);

if REV<FunSettings.PenaltyRev
    if REV<0
        REV = REV- 100*(2+abs((abs(REV)+1)/(FunSettings.PenaltyRev+0.1)));%10*abs(obj*(REV+5)); % use for revenue
    else
        REV = REV -10*(2+abs((10*REV+1)/(FunSettings.PenaltyRev+0.1)));%10*abs(obj*(REV+5)); % use for revenue
    end
end

if DAM>FunSettings.PenaltyDam
    %     DAM = DAM*(10 + abs(DAM));
end

obj(1)= -REV;
obj(2)= DAM;
% for i=1:length(DAMint)
%     obj(end+1)= DAMint(i);
% end

end

