function Output = Get_baseline_vals_Weibull(Data,surrogate,Prat,surrogate_type,WeibullIn,VBinSize,constTI)

% Getting the response for the case where the Power is constant to Prat.
% This modified version does it based on weibull distributions not on time
% series.
% The response includes basic loads, power and revenue as cumulative values that
% will be used later with the optimizer.
%
% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart


Ind.Shut = find(Data.Price<=0 | Data.V<4 | Data.V>24 | isnan(Data.Price) )';
Ind.On = setdiff(1:length(Data.V),Ind.Shut);
DataFilt.V = Data.V(Ind.On);
DataFilt.TI = Data.TI(Ind.On);
DataFilt.Price = Data.Price(Ind.On);

[~,Ewsp,~] = histcounts(DataFilt.V,4:VBinSize:24);
binsV = Ewsp(1:end-1)+diff(Ewsp)/2;
if isempty(WeibullIn)
    parmHat = wblfit(DataFilt.V); % 1 scale 2 shape
    BinPDF = wblpdf(binsV,parmHat(1),parmHat(2)) ;
else
    BinPDF = wblpdf(binsV,WeibullIn(1),WeibullIn(2)) ;
end

BinPDF = BinPDF/sum(BinPDF);  % normalize to 1
if surrogate_type == 1
    resp = Get_Points_From_Surrogate_noLoad(binsV',constTI*ones(length(binsV),1),Prat*ones(length(binsV),1),surrogate,'spline');
elseif surrogate_type == 2
    resp = Get_Points_From_Surrogate_GPR(binsV',constTI*ones(length(binsV),1),Prat*ones(length(binsV),1),surrogate);
end

% Weighted sum of all bins according to their normalized probability of occurence
Output.cum.Rev = sum(BinPDF'.*(resp.Power.mean/1000).*Data.Price(1));
Output.cum.Energy = sum(BinPDF'.*resp.Power.mean);
Output.cum.DAM.TBMx = sum(BinPDF'.*(resp.TBMx.mean.^4)*3600);
Output.cum.DAM.TBMy = sum(BinPDF'.*(resp.TBMy.mean.^4)*3600);
Output.cum.DAM.TBMz = sum(BinPDF'.*(resp.TBMz.mean.^4)*3600);
Output.cum.DAM.BRMx = sum(BinPDF'.*(resp.BRMx.mean.^10)*3600);
Output.cum.DAM.BRMy = sum(BinPDF'.*(resp.BRMy.mean.^10)*3600);
Output.cum.DAM.BRMz = sum(BinPDF'.*(resp.BRMz.mean.^10)*3600);
Output.cum.DAM.BROop = sum(BinPDF'.*(resp.BROop.mean.^10)*3600);
Output.cum.DAM.BRIp = sum(BinPDF'.*(resp.BRIp.mean.^10)*3600);
Output.cum.DAM.TTMx = sum(BinPDF'.*(resp.TTMx.mean.^4)*3600);
Output.cum.DAM.TTMy = sum(BinPDF'.*(resp.TTMy.mean.^4)*3600);
Output.cum.DAM.TTMz = sum(BinPDF'.*(resp.TTMz.mean.^4)*3600);
Output.cum.DAM.LSSMy = sum(BinPDF'.*(resp.LSSMy.mean.^4)*3600);
Output.cum.DAM.LSSTq = sum(BinPDF'.*(resp.LSSTq.mean.^4)*3600);




