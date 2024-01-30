function [Prat_out,optoptions,Thresholds] = Optimize_P_Shut_IBC_ConstantPrices_MultiSensor(VBinSize,WeibullIn,BaseLine,IBCflag,surrogate,TSData,surrogate_case,~,FunSettings)


% Optimize for revenue and damage accumulation for a constant price based on
% the general wsp distributio of the site. Output is the power level per wsp
% bin the revenue threshod to shut down and the IBC wsp/TI threshold to shut
% down
%
%
% - Create or load the weibull wsp distribution
% - Get the surrogate response based on the bin mean (or median) value
% - Objective function  evaluates Delta to damage and revenue compared to baseline for a given set of P levels per bin
% - Call the optimizer fmincon,ga,gamultiobj, multistart
%
% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart


for i =1:length(FunSettings.loadsens)
    if any(strcmp (FunSettings.loadsens{i},{'BROop';'BRIp';'BRMy';'BRMx';'BRMz'}))
        wohler(i) = 10;
    else
        wohler(i) = 4;
    end
    BaseDAM(i,1)  = BaseLine.cum.DAM.(FunSettings.loadsens{i});
end
BaseRev  = BaseLine.cum.Rev;


%% Get the wsp bins

% Identify the indices of the operating points
Ind.Shut = find(TSData.Price<=0 | TSData.V<4 | TSData.V>24 | isnan(TSData.Price) )';
Ind.On = setdiff(1:length(TSData.V),Ind.Shut);
Data.V = TSData.V(Ind.On);
Data.TI = TSData.TI(Ind.On);
Data.Price = TSData.Price(Ind.On);

% calculate edges bins and non fitted pdfs
[~,Ewsp,~] = histcounts(Data.V,4:VBinSize:24);
binsV = Ewsp(1:end-1)+diff(Ewsp)/2;
% [Nwsp,Ewsp,BinV] = histcounts(Data.V,4:VBinSize:24);
% curPDF =100*Nwsp/length(Data.V);

%% Calculate probability per bin of the current data
% get the pdf per bin assuming weibull fits
if isempty(WeibullIn)
    parmHat = wblfit(Data.V); % 1 scale 2 shape
    BinPDF = wblpdf(binsV,parmHat(1),parmHat(2)) ;
else
    BinPDF = wblpdf(binsV,WeibullIn(1),WeibullIn(2)) ;

end
BinPDF = BinPDF/sum(BinPDF);  % normalize to 100%

%% Assign bin bounds based on wsp
for iBin = 1:length(BinPDF)
    curV = binsV(iBin);
    if curV>=11
        lb(iBin) = 4.99;
        ub(iBin) = 11.01;
        x0(iBin) = 11;
    else
        lb(iBin) = 4.99;
        ub(iBin) = 10.01;
        x0(iBin) = 10;
    end
end

if FunSettings.constval.Shut==1
else
    x0(end+1) = FunSettings.constval.ShutVal;  % input for price thrweshold to shut down
    lb(end+1) = 0; % eur revenue threshold  (also a TI theshold here)
    ub(end+1) = 400; % eur revenue threshold  (also a TI theshold here)
end

if FunSettings.constval.IBC==1
else
    x0(end+1) = FunSettings.constval.IBCVal;  % input for IBC to start activating (it will not work with this sensor...)
    lb(end+1) = 12.99; %  m/s wsp to activate IBC (I could also add a TI threshold here)
    ub(end+1) = 24.01;  %  m/s wsp to activate IBC (I could also add a TI threshold here)
end


%% Optimization inputs
A = [];  % constrain the input according to actions for the region (A*x <= b)
b = [];  % constrain the actions according to actions for the region (A*x <= b)

Aeq = []; % linear equality constrains. Not applicable for me
beq = []; % linear equality constrains. Not applicable for me

options = optimoptions('fmincon', ...
    'DiffMinChange',0.5,...  % default 0
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
    'StepTolerance',1e-18, ... % default 1e-6 for all except interior point 1e-10
    'MaxFunctionEvaluations',15000, ...% default 3000
    'UseParallel',true...
    ...%     'ConstraintTolerance',1e-1,... % default 1e-6
    ...%     'OptimalityTolerance',1e-12...  % default 1e-9
    );

optionsga = optimoptions('ga',....
    'Display','iter', ...  % 'iter'   'off' 'diagnose'
    'UseParallel',true,...
    'FunctionTolerance', 1.0e-03,... %default 1e-6
    'MaxGenerations',100, ...  % default 100*length(inputs)
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
    'FunctionTolerance', 1.0e-06,... %default 1e-6
    'MaxGenerations',1500, ...  % default 100*length(inputs)
    'PopulationSize',100, ...   % defaut 200
    'CrossoverFraction',0.8, ... % default 0.8 how many pass to the next gen
    'MigrationFraction', 0.2, ... % default 0.2
    'PlotFcn',{@gaplotpareto} ...
    );

% Speed up
if surrogate_case==1
    [X1,Y1,Z1] = ndgrid(surrogate{1,1}.DataCnt.Dimensions.Dim1{2},surrogate{1, 1}.DataCnt.Dimensions.Dim2{2},surrogate{1,1}.DataCnt.Dimensions.Dim3{2});
end

if surrogate_case==1
    fun = @(x)func_obj(x,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.wDAM,FunSettings.wREV,FunSettings.loadsens,wohler,FunSettings.PenaltyRev,FunSettings.PenaltyDam,X1,Y1,Z1,Ind,Ewsp,IBCflag,FunSettings);
else
    fun = @(x)func_obj(x,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.wDAM,FunSettings.wREV,FunSettings.loadsens,wohler,FunSettings.PenaltyRev,FunSettings.PenaltyDam,[],[],[],Ind,Ewsp,IBCflag,FunSettings);
end

if surrogate_case==1
    funMulti = @(x)func_objMulti(x,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.loadsens,wohler,X1,Y1,Z1,Ind,Ewsp,IBCflag);
else
    funMulti = @(x)func_objMulti(x,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,FunSettings.loadsens,wohler,[],[],[],Ind,Ewsp,IBCflag);
end

%% Apply optimization

if strcmp(FunSettings.Method, 'fmincon')
    %%%     FMINCON       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    [x,fval,~,output,~]  = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,[],options);
    optoptions.Options = options;
    toc

elseif strcmp(FunSettings.Method, 'GlobalSearch')
    %%%   GLOBAL SEARCH   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    x = run(gs,problem);
    toc

elseif strcmp(FunSettings.Method, 'MultiStart')
    ms = MultiStart( ...
        'FunctionTolerance', 1e-3, ...  default 1e-6
        'Display', 'iter', ...
        'StartPointsToRun','all', ... default 'all'/ 'bounds'
        'UseParallel', true ...
        );
    problem = createOptimProblem('fmincon','x0',x0,'objective',fun,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',options);
    [x,fval,exitflag,output,solutions] = run(ms,problem,100);

elseif strcmp(FunSettings.Method, 'GA')
    %%% GENETIC ALGORITHM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    FunSettings.optimizerOptions = optionsga;
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
        x= xout(find(fval(:,1)==min(fval(allval,1))),:)'; %#ok<*FNDSB>
    else
        x= xout(find(fval(:,1)==min(fval(:,1))),:)';
    end
    toc
end


%% unpacking the bins to TS and assign outputs
x = round(x*10)/10;
for iT = 1:length(TSData.V)
    if TSData.Price(iT)<=0 || TSData.V(iT)<4 || TSData.V(iT)>24 || isnan(TSData.Price(iT))
        Prat_out(iT,1) =0;
    else
        row =  discretize(TSData.V(iT),Ewsp) ;
        Prat_out(iT,1) = x(row);
    end
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
optoptions.WeibullParametrs = WeibullIn;
optoptions.Vbins = binsV;


%% functions
function  obj = func_obj(x,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,wDAM,wRev,loadsens,wohler,PenaltyRev,PenaltyDam,X1,Y1,Z1,Ind,Ewsp,IBCflag,FunSettings) %#ok<INUSL>
% Calculate new cumulative values based on the weibull
% Compare to baseline
% Apply penatlies and weights

% separate the values to P and thresholds
if FunSettings.constval.Shut==0 && FunSettings.constval.IBC==0
    IBC_thres = x(end);
    Shut_thres = x(end-1);
elseif FunSettings.constval.IBC==1 && FunSettings.constval.Shut==1
    IBC_thres = FunSettings.constval.IBCVal;
    Shut_thres = FunSettings.constval.ShutVal;
    
elseif FunSettings.constval.IBC==0 && FunSettings.constval.Shut==1
    IBC_thres = x(end);
    Shut_thres = FunSettings.constval.ShutVal;
    x = x(1:end-1);
elseif FunSettings.constval.IBC==1 && FunSettings.constval.Shut==0
    IBC_thres = FunSettings.constval.IBCVal;
    Shut_thres = x(end);
    x = x(1:end-1);
end

% Assign IBC according to the new thresholds
Ind.newIBC = find(TSData.Price>0 & TSData.V>=IBC_thres & TSData.V<=24 & ~isnan(TSData.Price));
Ind.noIBC = setdiff(Ind.On,Ind.newIBC);
Ind.IBC   = setdiff(Ind.On,Ind.noIBC);

% Get P level per time block
P = zeros(length(TSData.V),1); % all power rating of TS
rows = discretize(TSData.V(Ind.On),Ewsp) ;
for i=1:length(Ind.On)
    P(Ind.On(i),1) = x(rows(i)); %#ok<*AGROW,*NASGU>
end

% Get the surrpgate response for each time step for both cases
SurOut.Power.mean =  zeros(length(TSData.V),1);
for i =1:length(loadsens)
    SurOut.(loadsens{i}).mean = zeros(length(TSData.V),1);
end
if surrogate_case == 1
    SurOut.Power.mean(Ind.noIBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.Power.mean,TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC),'spline'));
    for i =1:length(loadsens)
        % no IBC
        SurOut.(loadsens{i}).mean(Ind.noIBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,1}.DataCnt.(loadsens{i}).mean,TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC),'spline'));
        % IBC
        SurOut.(loadsens{i}).mean(Ind.IBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,2}.DataCnt.(loadsens{i}).mean,TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC),'spline'));
    end
    SurOut.Power.mean(Ind.IBC) = squeeze(interpn(X1,Y1,Z1,surrogate{1,2}.DataCnt.Power.mean,TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC),'spline'));
elseif surrogate_case == 2
    SurOut.Power.mean(Ind.noIBC) = predict(surrogate{1,1}.Output.Power,[TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC)]);
    for i =1:length(loadsens)
        % no IBC
        SurOut.(loadsens).mean(Ind.noIBC) = predict(surrogate{1,1}.Output.(loadsens),[TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC)]);
        % IBC
        SurOut.(loadsens).mean(Ind.IBC) = predict(surrogate{1,2}.Output.(loadsens),[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
    end
    SurOut.Power.mean(Ind.IBC) = predict(surrogate{1,2}.Output.Power,[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
end

% Find shut down bins according to the new threshold and remove them from the list to be evaluated
RevProj = (SurOut.Power.mean/1000).*TSData.Price; % projected revenue for the whole TS
Ind.ShutFinal = find(RevProj<Shut_thres);  % indices below threshold
SurOut.Power.mean(Ind.ShutFinal) = 0; % assign 0 power to these indices
for i =1:length(loadsens)  % assign 0 load to these indices
    SurOut.(loadsens{i}).mean(Ind.ShutFinal) = 0;
end

% Calculate new cummulative values
RevCum = sum(TSData.Price.*SurOut.Power.mean/1000);
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
if ~isempty(IndHigh)
    DAMint(IndHigh) = dlim; % to avoid load reductions dominating the single objective I used 10 for revenue and 20 for load
end
DAM  = mean(DAMint);   % here different types of weighting can be introduced

% Combining the objectives to one with weights
obj = -wRev*REV+wDAM*DAM;


if DAM > PenaltyDam
    %     obj = abs(obj*DAM); % use for revenue
    obj = obj+10; % use for load
elseif REV<-2*PenaltyRev
%     obj = 10*abs(obj*REV); % use for revenue
obj = obj+10*abs(1+abs(REV/(PenaltyRev+0.1))); % use for load
elseif REV<PenaltyRev
%     obj = abs(obj*REV); % use for revenue
    obj = obj+10*abs(1+abs(REV/(PenaltyRev+0.1)));  % use for load
end


function  obj = func_objMulti(x,surrogate,surrogate_case,BaseRev,BaseDAM,TSData,loadsens,wohler,X1,Y1,Z1,Ind,Ewsp,IBCflag)
% Calculate new cumulative values based on the weibull
% Compare to baseline
% Apply penatlies and weights

% separate the values to P and thresholds
if IBCflag==1
    IBC_thres  = 11;
else
    IBC_thres  = 25;
end

Shut_thres = 100;
% x = x(1:end-2);

% Assign IBC according to the new thresholds
Ind.newIBC = find(TSData.Price>0 & TSData.V>=IBC_thres & TSData.V<=24 & ~isnan(TSData.Price));
Ind.noIBC = setdiff(Ind.On,Ind.newIBC);
Ind.IBC   = setdiff(Ind.On,Ind.noIBC);

% Get P level per time block
P = zeros(length(TSData.V),1); % all power rating of TS
rows = discretize(TSData.V(Ind.On),Ewsp) ;
for i=1:length(Ind.On)
    P(Ind.On(i),1) = x(rows(i)); %#ok<*AGROW,*NASGU>
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
        SurOut.(loadsens).mean(Ind.noIBC) = predict(surrogate{1,1}.Output.(loadsens),[TSData.V(Ind.noIBC),TSData.TI(Ind.noIBC),P(Ind.noIBC)]);
        % IBC
        SurOut.(loadsens).mean(Ind.IBC) = predict(surrogate{1,2}.Output.(loadsens),[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
    end
    SurOut.Power.mean(Ind.IBC) = predict(surrogate{1,2}.Output.Power,[TSData.V(Ind.IBC),TSData.TI(Ind.IBC),P(Ind.IBC)]);
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
DAMint = 100*(DamCum./BaseDAM-1);
obj(1)= -REV;
for i=1:length(DAMint)
    obj(end+1)= DAMint(i);
end
