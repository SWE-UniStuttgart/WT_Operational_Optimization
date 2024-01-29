function Out = Get_Points_From_Surrogate_noLoad(Vreq,TIreq,Preq,sur_matrix,interp_type)
% close all

%% INPUT

% constTSR = load([surrogate_path '\constTSR' extension '.mat']);  % VxTIxP
% constTSR_IPC = load([surrogate_path '\constTSR_IPC' extension '.mat']); % VxTIxP
% lin70 = load([surrogate_path '\lin70' extension '.mat']);  % VxTIxP
% lin70_IPC = load([surrogate_path '\lin70_IPC' extension '.mat']);  % VxTIxP

% Vreq = 12; %Mean wind speed in m/s (DIM1) 4-24 range
% TIreq = 20;  % TI in %  (DIM2) 2-24 range
% Preq = 5:0.1:13;
%% Calculations

% Matlab native interpolation for gridded data linear,cubic,makima,spline

DataCnt = sur_matrix.DataCnt;
[X1,Y1,Z1] = ndgrid(DataCnt.Dimensions.Dim1{2},DataCnt.Dimensions.Dim2{2},DataCnt.Dimensions.Dim3{2});


% constTSR
Out.GenSpeedSTD.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.GenSpeedSTD.mean,Vreq,TIreq,Preq,interp_type));
Out.GenSpeedSTD.std   =squeeze(interpn(X1,Y1,Z1,DataCnt.GenSpeedSTD.std,Vreq,TIreq,Preq,interp_type));
Out.Power.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.Power.mean,Vreq,TIreq,Preq,interp_type));
Out.Power.std   = squeeze(interpn(X1,Y1,Z1,DataCnt.Power.std,Vreq,TIreq,Preq,interp_type));
Out.Energy.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.Energy.mean,Vreq,TIreq,Preq,interp_type));
Out.Energy.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.Energy.std,Vreq,TIreq,Preq,interp_type));
Out.GenTqSTD.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.GenTqSTD.mean,Vreq,TIreq,Preq,interp_type));
Out.GenTqSTD.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.GenTqSTD.std,Vreq,TIreq,Preq,interp_type));
Out.BlPitchSTD.mean   = squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchSTD.mean,Vreq,TIreq,Preq,interp_type));
Out.BlPitchSTD.std   = squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchSTD.std,Vreq,TIreq,Preq,interp_type));
Out.BlPitchTrav.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchTrav.mean,Vreq,TIreq,Preq,interp_type));
Out.BlPitchTrav.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchTrav.std,Vreq,TIreq,Preq,interp_type));
Out.TTDspFA.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTDspFA.mean,Vreq,TIreq,Preq,interp_type));
Out.TTDspFA.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTDspFA.std,Vreq,TIreq,Preq,interp_type));
Out.TBMx.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.TBMx.mean,Vreq,TIreq,Preq,interp_type));
Out.TBMx.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.TBMx.std,Vreq,TIreq,Preq,interp_type));
Out.TBMy.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.TBMy.mean,Vreq,TIreq,Preq,interp_type));
Out.TBMy.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.TBMy.std,Vreq,TIreq,Preq,interp_type));
Out.TBMz.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.TBMz.mean,Vreq,TIreq,Preq,interp_type));
Out.TBMz.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.TBMz.std,Vreq,TIreq,Preq,interp_type));
Out.BRMx.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.BRMx.mean,Vreq,TIreq,Preq,interp_type));
Out.BRMx.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.BRMx.std,Vreq,TIreq,Preq,interp_type));
Out.BRMy.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.BRMy.mean,Vreq,TIreq,Preq,interp_type));
Out.BRMy.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.BRMy.std,Vreq,TIreq,Preq,interp_type));
Out.BRMz.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.BRMz.mean,Vreq,TIreq,Preq,interp_type));
Out.BRMz.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.BRMz.std,Vreq,TIreq,Preq,interp_type));
Out.BROop.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.BROop.mean,Vreq,TIreq,Preq,interp_type));
Out.BROop.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.BROop.std,Vreq,TIreq,Preq,interp_type));
Out.BRIp.mean   = squeeze(interpn(X1,Y1,Z1,DataCnt.BRIp.mean,Vreq,TIreq,Preq,interp_type));
Out.BRIp.std   = squeeze(interpn(X1,Y1,Z1,DataCnt.BRIp.std,Vreq,TIreq,Preq,interp_type));
Out.TTMx.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTMx.mean,Vreq,TIreq,Preq,interp_type));
Out.TTMx.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTMx.std,Vreq,TIreq,Preq,interp_type));
Out.TTMy.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTMy.mean,Vreq,TIreq,Preq,interp_type));
Out.TTMy.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTMy.std,Vreq,TIreq,Preq,interp_type));
Out.TTMz.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTMz.mean,Vreq,TIreq,Preq,interp_type));
Out.TTMz.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.TTMz.std,Vreq,TIreq,Preq,interp_type));
Out.LSSMy.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMy.mean,Vreq,TIreq,Preq,interp_type));
Out.LSSMy.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMy.std,Vreq,TIreq,Preq,interp_type));
Out.LSSMz.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMz.mean,Vreq,TIreq,Preq,interp_type));
Out.LSSMz.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMz.std,Vreq,TIreq,Preq,interp_type));
Out.LSSTq.mean  = squeeze(interpn(X1,Y1,Z1,DataCnt.LSSTq.mean,Vreq,TIreq,Preq,interp_type));
Out.LSSTq.std  = squeeze(interpn(X1,Y1,Z1,DataCnt.LSSTq.std,Vreq,TIreq,Preq,interp_type));

% Out.GenSpeedSTD(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.GenSpeedSTD.mean,Vreq,TIreq,Preq,interp_type));
% Out.GenSpeedSTD(:,2) =squeeze(interpn(X1,Y1,Z1,DataCnt.GenSpeedSTD.std,Vreq,TIreq,Preq,interp_type));
% Out.Power(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.Power.mean,Vreq,TIreq,Preq,interp_type));
% Out.Power(:,2) = squeeze(interpn(X1,Y1,Z1,DataCnt.Power.std,Vreq,TIreq,Preq,interp_type));
% Out.Energy(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.Energy.mean,Vreq,TIreq,Preq,interp_type));
% Out.Energy(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.Energy.std,Vreq,TIreq,Preq,interp_type));
% Out.GenTqSTD(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.GenTqSTD.mean,Vreq,TIreq,Preq,interp_type));
% Out.GenTqSTD(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.GenTqSTD.std,Vreq,TIreq,Preq,interp_type));
% Out.BlPitchSTD(:,1) = squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchSTD.mean,Vreq,TIreq,Preq,interp_type));
% Out.BlPitchSTD(:,2) = squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchSTD.std,Vreq,TIreq,Preq,interp_type));
% Out.BlPitchTrav(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchTrav.mean,Vreq,TIreq,Preq,interp_type));
% Out.BlPitchTrav(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.BlPitchTrav.std,Vreq,TIreq,Preq,interp_type));
% Out.TTDspFA(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTDspFA.mean,Vreq,TIreq,Preq,interp_type));
% Out.TTDspFA(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTDspFA.std,Vreq,TIreq,Preq,interp_type));
% Out.TBMx(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.TBMx.mean,Vreq,TIreq,Preq,interp_type));
% Out.TBMx(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.TBMx.std,Vreq,TIreq,Preq,interp_type));
% Out.TBMy(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.TBMy.mean,Vreq,TIreq,Preq,interp_type));
% Out.TBMy(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.TBMy.std,Vreq,TIreq,Preq,interp_type));
% Out.TBMz(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.TBMz.mean,Vreq,TIreq,Preq,interp_type));
% Out.TBMz(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.TBMz.std,Vreq,TIreq,Preq,interp_type));
% Out.BRMx(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.BRMx.mean,Vreq,TIreq,Preq,interp_type));
% Out.BRMx(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.BRMx.std,Vreq,TIreq,Preq,interp_type));
% Out.BRMy(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.BRMy.mean,Vreq,TIreq,Preq,interp_type));
% Out.BRMy(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.BRMy.std,Vreq,TIreq,Preq,interp_type));
% Out.BRMz(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.BRMz.mean,Vreq,TIreq,Preq,interp_type));
% Out.BRMz(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.BRMz.std,Vreq,TIreq,Preq,interp_type));
% Out.BROop(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.BROop.mean,Vreq,TIreq,Preq,interp_type));
% Out.BROop(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.BROop.std,Vreq,TIreq,Preq,interp_type));
% Out.BRIp(:,1) = squeeze(interpn(X1,Y1,Z1,DataCnt.BRIp.mean,Vreq,TIreq,Preq,interp_type));
% Out.BRIp(:,2) = squeeze(interpn(X1,Y1,Z1,DataCnt.BRIp.std,Vreq,TIreq,Preq,interp_type));
% Out.TTMx(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTMx.mean,Vreq,TIreq,Preq,interp_type));
% Out.TTMx(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTMx.std,Vreq,TIreq,Preq,interp_type));
% Out.TTMy(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTMy.mean,Vreq,TIreq,Preq,interp_type));
% Out.TTMy(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTMy.std,Vreq,TIreq,Preq,interp_type));
% Out.TTMz(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTMz.mean,Vreq,TIreq,Preq,interp_type));
% Out.TTMz(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.TTMz.std,Vreq,TIreq,Preq,interp_type));
% Out.LSSMy(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMy.mean,Vreq,TIreq,Preq,interp_type));
% Out.LSSMy(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMy.std,Vreq,TIreq,Preq,interp_type));
% Out.LSSMz(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMz.mean,Vreq,TIreq,Preq,interp_type));
% Out.LSSMz(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.LSSMz.std,Vreq,TIreq,Preq,interp_type));
% Out.LSSTq(:,1)= squeeze(interpn(X1,Y1,Z1,DataCnt.LSSTq.mean,Vreq,TIreq,Preq,interp_type));
% Out.LSSTq(:,2)= squeeze(interpn(X1,Y1,Z1,DataCnt.LSSTq.std,Vreq,TIreq,Preq,interp_type));


end