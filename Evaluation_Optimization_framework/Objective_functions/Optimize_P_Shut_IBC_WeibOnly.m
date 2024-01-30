function [Prat_out,optoptions] = Optimize_P_Shut_IBC_WeibOnly(TSData,FunSettings)


% Getting an omld results based on constant prices with weibull
% optimization and applying it on a fluctuating prices
%
% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart

%% Get the wsp bins

% Identify the indices of the operating points
Ind.Shut = find(TSData.Price<=0 | TSData.V<4 | TSData.V>24 | isnan(TSData.Price) )';
Ind.On = setdiff(1:length(TSData.V),Ind.Shut);
Data.V = TSData.V(Ind.On);
Data.TI = TSData.TI(Ind.On);
Data.Price = TSData.Price(Ind.On);

OldRes = load(FunSettings.ResFile);
VBinSize = OldRes.Output.options.control.logic.VbinSize;
% calculate bin edges and centers 
[~,Ewsp,~] = histcounts(Data.V,4:VBinSize:24);
binsV = Ewsp(1:end-1)+diff(Ewsp)/2;

%% unpacking the bins to TS and assign outputs
x = OldRes.Output.options.control.FunSettings.optimizer{1, 1}.Values(1:end-2);
x = round(x*10)/10;
for iT = 1:length(TSData.V)
    if TSData.Price(iT)<=0 || TSData.V(iT)<4 || TSData.V(iT)>24 || isnan(TSData.Price(iT))
        Prat_out(iT,1) =0; 
    else
        row =  discretize(TSData.V(iT),Ewsp) ;
        Prat_out(iT,1) = x(row); %#ok<*AGROW> 
    end
end

Thresholds.IBC.WSP = OldRes.Output.options.control.FunSettings.optimizer{1, 1}.Thresholds.IBC.WSP  ;
Thresholds.Shut.Rev = OldRes.Output.options.control.FunSettings.optimizer{1, 1}.Thresholds.Shut.Rev  ;
optoptions.Values = x;
optoptions.Thresholds = Thresholds;
optoptions.Vbins = binsV;


