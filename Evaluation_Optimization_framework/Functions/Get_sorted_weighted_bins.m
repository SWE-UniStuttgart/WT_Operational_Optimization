function Out = Get_sorted_weighted_bins(Baseline,Bins)

% sorting the bins according to their contribution to revenue and DEL

% Go through each bin and addthe value of damage and revenue
for iV = 1:size(Bins.Prdisc,1)
    for iPr = 1:size(Bins.Prdisc,2)
        if sum(Bins.Pbin{iV,iPr}(:,1)) ==0 %find empty bins
            RevBase(iV,iPr) = 0; %#ok<*AGROW>
            DAMBase.TBMx(iV,iPr)  = 0;
            DAMBase.TBMy(iV,iPr)  = 0;
            DAMBase.TBMz(iV,iPr)  = 0;
            DAMBase.BRMx(iV,iPr)  = 0;
            DAMBase.BRMy(iV,iPr)  = 0;
            DAMBase.BRMz(iV,iPr)  = 0;
            DAMBase.BROop(iV,iPr) = 0;
            DAMBase.BRIp(iV,iPr) = 0;
            DAMBase.TTMy(iV,iPr)  = 0;
            DAMBase.LSSMy(iV,iPr) = 0;
            DAMBase.LSSTq(iV,iPr) = 0;
        else
            % go through the indices (indicating time steps of the baseline TS) included in the bins and add their values to the bin
            % This does the weighting automatically
            RevBase(iV,iPr) =  sum(Baseline.Rev(Bins.Pbin{iV,iPr}(:,1)));
            DAMBase.TBMx(iV,iPr)  = sum( Baseline.DAM.TBMx(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.TBMy(iV,iPr)  = sum( Baseline.DAM.TBMy(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.TBMz(iV,iPr)  = sum( Baseline.DAM.TBMz(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.BRMx(iV,iPr)  = sum( Baseline.DAM.BRMx(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.BRMy(iV,iPr)  = sum( Baseline.DAM.BRMy(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.BRMz(iV,iPr)  = sum( Baseline.DAM.BRMz(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.BROop(iV,iPr) = sum( Baseline.DAM.BROop(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.BRIp(iV,iPr) = sum( Baseline.DAM.BRIp(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.TTMy(iV,iPr)  = sum( Baseline.DAM.TTMy(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.LSSMy(iV,iPr) = sum( Baseline.DAM.LSSMy(Bins.Pbin{iV,iPr}(:,1)) );
            DAMBase.LSSTq(iV,iPr) = sum( Baseline.DAM.LSSTq(Bins.Pbin{iV,iPr}(:,1)) );
        end
    end
end
Out.Rev_Base_ND = RevBase/sum(sum(RevBase));
Out.Rev_Base = RevBase;
Out.DAMBase  = DAMBase;
Out.DAMBase_ND.TBMx = DAMBase.TBMx/sum(sum(DAMBase.TBMx));
Out.DAMBase_ND.TBMy = DAMBase.TBMy/sum(sum(DAMBase.TBMy));
Out.DAMBase_ND.TBMz = DAMBase.TBMz/sum(sum(DAMBase.TBMz));
Out.DAMBase_ND.BRMx = DAMBase.BRMx/sum(sum(DAMBase.BRMx));
Out.DAMBase_ND.BRMy = DAMBase.BRMy/sum(sum(DAMBase.BRMy));
Out.DAMBase_ND.BRMz = DAMBase.BRMz/sum(sum(DAMBase.BRMz));
Out.DAMBase_ND.BROop = DAMBase.BROop/sum(sum(DAMBase.BROop));
Out.DAMBase_ND.BRIp = DAMBase.BROop/sum(sum(DAMBase.BRIp));
Out.DAMBase_ND.TTMy  = DAMBase.TTMy/sum(sum(DAMBase.TTMy));
Out.DAMBase_ND.LSSMy = DAMBase.LSSMy/sum(sum(DAMBase.LSSMy));
Out.DAMBase_ND.LSSTq = DAMBase.LSSTq/sum(sum(DAMBase.LSSTq));

[Out.Rev_Base_sort(:,1),Out.Rev_Base_sort(:,2)] = sort(Out.Rev_Base(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.TBMx_sort(:,1),Out.DAMBase.TBMx_sort(:,2)] = sort(Out.DAMBase.TBMx(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.TBMy_sort(:,1),Out.DAMBase.TBMy_sort(:,2)] = sort(Out.DAMBase.TBMy(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.TBMz_sort(:,1),Out.DAMBase.TBMz_sort(:,2)] = sort(Out.DAMBase.TBMz(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.BRMx_sort(:,1),Out.DAMBase.BRMx_sort(:,2)] = sort(Out.DAMBase.BRMx(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.BRMy_sort(:,1),Out.DAMBase.BRMy_sort(:,2)] = sort(Out.DAMBase.BRMy(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.BRMz_sort(:,1),Out.DAMBase.BRMz_sort(:,2)] = sort(Out.DAMBase.BRMz(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.BROop_sort(:,1),Out.DAMBase.BROop_sort(:,2)] = sort(Out.DAMBase.BROop(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.BRIp_sort(:,1),Out.DAMBase.BRIp_sort(:,2)] = sort(Out.DAMBase.BRIp(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.TTMy_sort(:,1),Out.DAMBase.TTMy_sort(:,2)]   = sort(Out.DAMBase.TTMy(:),'descend'); % 1 val 2 linear index
[Out.DAMBase.LSSMy_sort(:,1),Out.DAMBase.LSSMy_sort(:,2)] = sort(Out.DAMBase.LSSMy(:),'descend'); % 1 val 2 linearindex
[Out.DAMBase.LSSTq_sort(:,1),Out.DAMBase.LSSTq_sort(:,2)] = sort(Out.DAMBase.LSSTq(:),'descend'); % 1 val 2 linearindex






