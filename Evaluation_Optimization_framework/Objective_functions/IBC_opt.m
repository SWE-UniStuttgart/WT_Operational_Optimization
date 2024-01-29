% The optimization function for IBC only. It takes the current
% conditions and evaluates whether to switch IBC on.



function [IBCon] = IBC_opt(V,Power,TI,FunSettings)

if TI<=FunSettings.TI_Thres_IBC
    if V>FunSettings.VThresLowTI_IBC
        IBCon = 1;
    else
        IBCon = 0;
    end
else
    if V>FunSettings.VThresHighTI_IBC
        IBCon = 1;
    else
        IBCon = 0;
    end
end



