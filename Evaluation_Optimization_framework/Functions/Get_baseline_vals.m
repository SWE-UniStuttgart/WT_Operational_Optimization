function Output = Get_baseline_vals(Data,surrogate,Prat,surrogate_type)

% Getting the response for the case where the Power is constant to Prat.
% The response includes basic loads, power and revenue as time series that
% will be used later by the optimizer.
%
% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart


% Pre-calculate GPR response it is much faster

IndShut = find(Data.V <4 | Data.V >24 | isnan(Data.Price) | Data.Price<=0);
IndOn = setdiff([1:length(Data.V)]',IndShut);

Operating.V = Data.V(IndOn);
Operating.TI = Data.TI(IndOn);
Operating.Price = Data.Price(IndOn);

if surrogate_type == 1
    resp = Get_Points_From_Surrogate_noLoad(Operating.V,Operating.TI,Prat*ones(length(Operating.V),1),surrogate,'spline');
elseif surrogate_type == 2
    resp = Get_Points_From_Surrogate_GPR(Operating.V,Operating.TI,Prat*ones(length(Operating.V),1),surrogate);
end


Output.Rev(IndOn) = resp.Power.mean/1000.*Operating.Price;
Output.Energy(IndOn) = resp.Power.mean;
Output.DAM.TBMx(IndOn)  = (resp.TBMx.mean.^4)*3600;
Output.DAM.TBMy(IndOn)  = (resp.TBMy.mean.^4)*3600;
Output.DAM.TBMz(IndOn)  = (resp.TBMz.mean.^4)*3600;
Output.DAM.BRMx(IndOn)  = (resp.BRMx.mean.^10)*3600;
Output.DAM.BRMy(IndOn)  = (resp.BRMy.mean.^10)*3600;
Output.DAM.BRMz(IndOn)  = (resp.BRMz.mean.^10)*3600;
Output.DAM.BROop(IndOn) = (resp.BROop.mean.^10)*3600;
Output.DAM.BRIp(IndOn) = (resp.BRIp.mean.^10)*3600;
Output.DAM.TTMx(IndOn)  = (resp.TTMx.mean.^4)*3600;
Output.DAM.TTMy(IndOn)  = (resp.TTMy.mean.^4)*3600;
Output.DAM.TTMz(IndOn)  = (resp.TTMz.mean.^4)*3600;
Output.DAM.LSSMy(IndOn) = (resp.LSSMy.mean.^4)*3600;
Output.DAM.LSSTq(IndOn) = (resp.LSSTq.mean.^4)*3600;

Output.Rev(IndShut) = 0;
Output.Energy(IndShut) = 0;
Output.DAM.TBMx(IndShut)  = 0;
Output.DAM.TBMy(IndShut)  = 0;
Output.DAM.TBMz(IndShut)  = 0;
Output.DAM.BRMx(IndShut)  = 0;
Output.DAM.BRMy(IndShut)  = 0;
Output.DAM.BRMz(IndShut)  = 0;
Output.DAM.BROop(IndShut) = 0;
Output.DAM.BRIp(IndShut) = 0;
Output.DAM.TTMx(IndShut)  = 0;
Output.DAM.TTMy(IndShut)  = 0;
Output.DAM.TTMz(IndShut)  = 0;
Output.DAM.LSSMy(IndShut) = 0;
Output.DAM.LSSTq(IndShut) = 0;

Output.cum.Rev = sum(resp.Power.mean/1000.*Operating.Price);
Output.cum.Energy= sum(resp.Power.mean);
Output.cum.DAM.TBMx  = sum(Output.DAM.TBMx);
Output.cum.DAM.TBMy  = sum(Output.DAM.TBMy);
Output.cum.DAM.TBMz  = sum(Output.DAM.TBMz);
Output.cum.DAM.BRMx  = sum(Output.DAM.BRMx);
Output.cum.DAM.BRMy  = sum(Output.DAM.BRMy);
Output.cum.DAM.BRMz  = sum(Output.DAM.BRMz);
Output.cum.DAM.BROop = sum(Output.DAM.BROop);
Output.cum.DAM.BRIp = sum(Output.DAM.BRIp);
Output.cum.DAM.TTMx  = sum(Output.DAM.TTMx);
Output.cum.DAM.TTMy  = sum(Output.DAM.TTMy);
Output.cum.DAM.TTMz  = sum(Output.DAM.TTMz);
Output.cum.DAM.LSSMy = sum(Output.DAM.LSSMy);
Output.cum.DAM.LSSTq = sum(Output.DAM.LSSTq);




