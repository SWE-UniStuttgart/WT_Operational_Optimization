% The optimization function for shutting down only. It takes the current
% conditions and evaluates whether to shut down.

% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart

function [Shut]=Shut_down_opt(V,Pr,TI,Pcurve,FunSettings)


cur.P = interp1(Pcurve.power(:,1),Pcurve.power(:,2),V);
cur.L = interp1(Pcurve.load(:,1),Pcurve.load(:,2),V);
cur.Rev = cur.P*Pr;

if TI<=FunSettings.TI_Thres_Shut
    if cur.Rev<=FunSettings.revThresLowTI_Shut
        Shut = 1;
    else
        Shut = 0;
    end
else
    if cur.Rev<=FunSettings.revThresHighTI_Shut
        Shut = 1;
    else
        Shut = 0;
    end
end






