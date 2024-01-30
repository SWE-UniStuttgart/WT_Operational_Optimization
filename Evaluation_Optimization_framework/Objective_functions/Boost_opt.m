% The optimization function for power boosting. It is used only in arated
% annd above wind speeds to determine whether to boost and at what level.
% There is no optimization involved it is a static feeforward with thresholds on
% conditions
%
% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart

function [Plevel,FunSettings] = Boost_opt(V,Pr,TI,Pcurve,FunSettings)

cur.P = interp1(Pcurve.power(:,1),Pcurve.power(:,2),V);
cur.L = interp1(Pcurve.load(:,1),Pcurve.load(:,2),V);
cur.Rev = cur.P*Pr;

if cur.Rev>= FunSettings.revThresHigh2_Boost && V>FunSettings.VThresHigh2_Boost
    Plevel = FunSettings.Pmax2;
elseif cur.Rev>= FunSettings.revThresHigh_Boost && V>FunSettings.VThresHigh_Boost
     Plevel = FunSettings.Pmax;
elseif cur.Rev>= FunSettings.revThresHigh_Boost && V<FunSettings.VThresHigh_Boost
     Plevel = FunSettings.Ptrans;
elseif cur.Rev> FunSettings.revThresMid_Boost && V>FunSettings.VThresHigh_Boost
    Plevel = FunSettings.Pmid;
elseif cur.Rev<= FunSettings.revThresHigh_Boost && V<FunSettings.VThresHigh_Boost
    Plevel = FunSettings.Plow;
else
    if TI< FunSettings.TI_Thres_Boost %#ok<IFBDUP> 
        Plevel = 10;
    else
%         Plevel = 10+ (cur.Rev-FunSettings.revThresMid_Boost)*(13-10)/(FunSettings.revThresHigh_Boost-FunSettings.revThresMid_Boost) ;%ramp
        Plevel = 10; %ramp
    end

end


%% old implementation
% cur.P = interp1(Pcurve.power(:,1),Pcurve.power(:,2),V);
% cur.L = interp1(Pcurve.load(:,1),Pcurve.load(:,2),V);
% cur.Rev = cur.P*Pr;
% 
% if cur.Rev>= FunSettings.revThresHigh_Boost
%      Plevel = 13;
% elseif cur.Rev< FunSettings.revThresMid_Boost
%     Plevel = 10;
% else
%     if TI< FunSettings.TI_Thres_Boost
%         Plevel = 13;
%     else
%         Plevel = 10+ (cur.Rev-FunSettings.revThresMid_Boost)*(13-10)/(FunSettings.revThresHigh_Boost-FunSettings.revThresMid_Boost) ;%ramp
%     end
% 
% end



