clc
clearvars
close all


% Script that takes reults and combines them to one to see the total
% effect. The files Areadded in the order defined

%% INputs

pathbase= 'YourPath\';
pathopt ='YourPath\';


OptFiles = {...
    'File1'
    'File2'
    'File3'
    }';


BaseFile = {
    'Baseline_DE_all_spl_const80.mat'
    };



%% Calculations

% Load the files
for i = 1:length(OptFiles)
    Opt{i,1} = load([pathopt OptFiles{i}]);
end
Baseline = load([pathbase BaseFile{1}]);

Output.Time{1,2} = [];

Output.V                  =[];
Output.Price              =[];
Output.TI                 =[];

% Add the inst and cum values
Output.inst.Cntr          = {};
Output.inst.Prat          = [];
Output.inst.IBC           = [];
Output.inst.Energy        =[];
Output.inst.Revenue       =[];
Output.inst.BlPitchTrav   =[];
Output.inst.BlPitchSTD    =[];
Output.inst.GenTqSTD      =[];
Output.inst.GenSpeedSTD   =  [];

Output.inst.DEL.TBMx      =[];
Output.inst.DEL.TBMy      =[];
Output.inst.DEL.TBMz      =[];
Output.inst.DEL.BRMx      =[];
Output.inst.DEL.BRMy      =[];
Output.inst.DEL.BRMz      =[];
Output.inst.DEL.BROop      =[];
Output.inst.DEL.BRIp      =[];
Output.inst.DEL.TTMx      =[];
Output.inst.DEL.TTMy      =[];
Output.inst.DEL.TTMz      =[];
Output.inst.DEL.LSSMy      =[];
Output.inst.DEL.LSSMz      =[];
Output.inst.DEL.LSSTq      =[];

Output.inst.DAM.TBMx      =[];
Output.inst.DAM.TBMy      =[];
Output.inst.DAM.TBMz      =[];
Output.inst.DAM.BRMx      =[];
Output.inst.DAM.BRMy      =[];
Output.inst.DAM.BRMz      =[];
Output.inst.DAM.BROop      =[];
Output.inst.DAM.BRIp      =[];
Output.inst.DAM.TTMx      =[];
Output.inst.DAM.TTMy      =[];
Output.inst.DAM.TTMz      =[];
Output.inst.DAM.LSSMy      =[];
Output.inst.DAM.LSSMz      =[];
Output.inst.DAM.LSSTq      =[];
% Output.inst.             =[];

Output.cum.Energy        = [];
Output.cum.Revenue       = [];
Output.cum.BlPitchTrav   = [];

Output.cum.DAM.TBMx      = [];
Output.cum.DAM.TBMy      = [];
Output.cum.DAM.TBMz      = [];
Output.cum.DAM.BRMx      = [];
Output.cum.DAM.BRMy      = [];
Output.cum.DAM.BRMz      = [];
Output.cum.DAM.BROop     = [];
Output.cum.DAM.BRIp      = [];
Output.cum.DAM.TTMx      = [];
Output.cum.DAM.TTMy      = [];
Output.cum.DAM.TTMz      = [];
Output.cum.DAM.LSSMy     = [];
Output.cum.DAM.LSSMz     = [];
Output.cum.DAM.LSSTq     = [];

for i = 1:length(OptFiles)
    if i ==1
        Output.Time{1,1}               = [Opt{i}.Output.Time{1,1}];
    else
        Output.Time{1,1}               = [Output.Time{1,1}             ; Opt{i}.Output.Time{1,1}];
    end
    Output.Time{1,2}               = [Output.Time{1,2}             ; Opt{i}.Output.Time{1,2}];


    Output.V                  = [Output.V                ; Opt{i}.Output.V ];
    Output.Price              = [Output.Price            ; Opt{i}.Output.Price];
    Output.TI                 = [Output.TI               ; Opt{i}.Output.TI];

    Output.inst.Cntr          = [Output.inst.Cntr        ; Opt{i}.Output.inst.Cntr           ];
    Output.inst.Prat          = [Output.inst.Prat        ; Opt{i}.Output.inst.Prat           ];
    Output.inst.IBC           = [Output.inst.IBC         ; Opt{i}.Output.inst.IBC            ];
    Output.inst.Energy        = [Output.inst.Energy      ; Opt{i}.Output.inst.Energy         ];
    Output.inst.Revenue       = [Output.inst.Revenue     ; Opt{i}.Output.inst.Revenue        ];
    Output.inst.BlPitchTrav   = [Output.inst.BlPitchTrav ; Opt{i}.Output.inst.BlPitchTrav    ];
    Output.inst.BlPitchSTD    = [Output.inst.BlPitchSTD  ; Opt{i}.Output.inst.BlPitchSTD     ];
    Output.inst.GenTqSTD      = [Output.inst.GenTqSTD    ; Opt{i}.Output.inst.GenTqSTD       ];
    Output.inst.GenSpeedSTD   = [Output.inst.GenSpeedSTD ; Opt{i}.Output.inst.GenSpeedSTD    ];

    Output.inst.DEL.TBMx      = [Output.inst.DEL.TBMx    ; Opt{i}.Output.inst.DEL.TBMx       ];
    Output.inst.DEL.TBMy      = [Output.inst.DEL.TBMy    ; Opt{i}.Output.inst.DEL.TBMy       ];
    Output.inst.DEL.TBMz      = [Output.inst.DEL.TBMz    ; Opt{i}.Output.inst.DEL.TBMz       ];
    Output.inst.DEL.BRMx      = [Output.inst.DEL.BRMx    ; Opt{i}.Output.inst.DEL.BRMx       ];
    Output.inst.DEL.BRMy      = [Output.inst.DEL.BRMy    ; Opt{i}.Output.inst.DEL.BRMy       ];
    Output.inst.DEL.BRMz      = [Output.inst.DEL.BRMz    ; Opt{i}.Output.inst.DEL.BRMz       ];
    Output.inst.DEL.BROop     = [Output.inst.DEL.BROop   ; Opt{i}.Output.inst.DEL.BROop      ];
    Output.inst.DEL.BRIp      = [Output.inst.DEL.BRIp    ; Opt{i}.Output.inst.DEL.BRIp       ];
    Output.inst.DEL.TTMx      = [Output.inst.DEL.TTMx    ; Opt{i}.Output.inst.DEL.TTMx       ];
    Output.inst.DEL.TTMy      = [Output.inst.DEL.TTMy    ; Opt{i}.Output.inst.DEL.TTMy       ];
    Output.inst.DEL.TTMz      = [Output.inst.DEL.TTMz    ; Opt{i}.Output.inst.DEL.TTMz       ];
    Output.inst.DEL.LSSMy     = [Output.inst.DEL.LSSMy   ; Opt{i}.Output.inst.DEL.LSSMy      ];
    Output.inst.DEL.LSSMz     = [Output.inst.DEL.LSSMz   ; Opt{i}.Output.inst.DEL.LSSMz      ];
    Output.inst.DEL.LSSTq     = [Output.inst.DEL.LSSTq   ; Opt{i}.Output.inst.DEL.LSSTq      ];

    Output.inst.DAM.TBMx      = [Output.inst.DAM.TBMx    ; Opt{i}.Output.inst.DAM.TBMx       ];
    Output.inst.DAM.TBMy      = [Output.inst.DAM.TBMy    ; Opt{i}.Output.inst.DAM.TBMy       ];
    Output.inst.DAM.TBMz      = [Output.inst.DAM.TBMz    ; Opt{i}.Output.inst.DAM.TBMz       ];
    Output.inst.DAM.BRMx      = [Output.inst.DAM.BRMx    ; Opt{i}.Output.inst.DAM.BRMx       ];
    Output.inst.DAM.BRMy      = [Output.inst.DAM.BRMy    ; Opt{i}.Output.inst.DAM.BRMy       ];
    Output.inst.DAM.BRMz      = [Output.inst.DAM.BRMz    ; Opt{i}.Output.inst.DAM.BRMz       ];
    Output.inst.DAM.BROop     = [Output.inst.DAM.BROop   ; Opt{i}.Output.inst.DAM.BROop      ];
    Output.inst.DAM.BRIp      = [Output.inst.DAM.BRIp    ; Opt{i}.Output.inst.DAM.BRIp       ];
    Output.inst.DAM.TTMx      = [Output.inst.DAM.TTMx    ; Opt{i}.Output.inst.DAM.TTMx       ];
    Output.inst.DAM.TTMy      = [Output.inst.DAM.TTMy    ; Opt{i}.Output.inst.DAM.TTMy       ];
    Output.inst.DAM.TTMz      = [Output.inst.DAM.TTMz    ; Opt{i}.Output.inst.DAM.TTMz       ];
    Output.inst.DAM.LSSMy     = [Output.inst.DAM.LSSMy   ; Opt{i}.Output.inst.DAM.LSSMy      ];
    Output.inst.DAM.LSSMz     = [Output.inst.DAM.LSSMz   ; Opt{i}.Output.inst.DAM.LSSMz      ];
    Output.inst.DAM.LSSTq     = [Output.inst.DAM.LSSTq   ; Opt{i}.Output.inst.DAM.LSSTq      ];

    if i==1
        Output.cum.Energy        = [Output.cum.Energy ; Opt{i}.Output.cum.Energy];
        Output.cum.Revenue       = [Output.cum.Revenue ;  Opt{i}.Output.cum.Revenue];
        Output.cum.BlPitchTrav   = [Output.cum.BlPitchTrav ; Opt{i}.Output.cum.BlPitchTrav];

        Output.cum.DAM.TBMx      = [Output.cum.DAM.TBMx    ; Opt{i}.Output.cum.DAM.TBMx       ];
        Output.cum.DAM.TBMy      = [Output.cum.DAM.TBMy    ; Opt{i}.Output.cum.DAM.TBMy       ];
        Output.cum.DAM.TBMz      = [Output.cum.DAM.TBMz    ; Opt{i}.Output.cum.DAM.TBMz       ];
        Output.cum.DAM.BRMx      = [Output.cum.DAM.BRMx    ; Opt{i}.Output.cum.DAM.BRMx       ];
        Output.cum.DAM.BRMy      = [Output.cum.DAM.BRMy    ; Opt{i}.Output.cum.DAM.BRMy       ];
        Output.cum.DAM.BRMz      = [Output.cum.DAM.BRMz    ; Opt{i}.Output.cum.DAM.BRMz       ];
        Output.cum.DAM.BROop     = [Output.cum.DAM.BROop   ; Opt{i}.Output.cum.DAM.BROop      ];
        Output.cum.DAM.BRIp      = [Output.cum.DAM.BRIp    ; Opt{i}.Output.cum.DAM.BRIp       ];
        Output.cum.DAM.TTMx      = [Output.cum.DAM.TTMx    ; Opt{i}.Output.cum.DAM.TTMx       ];
        Output.cum.DAM.TTMy      = [Output.cum.DAM.TTMy    ; Opt{i}.Output.cum.DAM.TTMy       ];
        Output.cum.DAM.TTMz      = [Output.cum.DAM.TTMz    ; Opt{i}.Output.cum.DAM.TTMz       ];
        Output.cum.DAM.LSSMy     = [Output.cum.DAM.LSSMy   ; Opt{i}.Output.cum.DAM.LSSMy      ];
        Output.cum.DAM.LSSMz     = [Output.cum.DAM.LSSMz   ; Opt{i}.Output.cum.DAM.LSSMz      ];
        Output.cum.DAM.LSSTq     = [Output.cum.DAM.LSSTq   ; Opt{i}.Output.cum.DAM.LSSTq      ];
    else
        Output.cum.Energy        = [Output.cum.Energy      ; Output.cum.Energy(end)    +  Opt{i}.Output.cum.Energy         ];
        Output.cum.Revenue       = [Output.cum.Revenue     ;Output.cum.Revenue(end) +   Opt{i}.Output.cum.Revenue      ];
        Output.cum.BlPitchTrav   = [Output.cum.BlPitchTrav ;Output.cum.BlPitchTrav(end) +  Opt{i}.Output.cum.BlPitchTrav    ];

        Output.cum.DAM.TBMx      = [Output.cum.DAM.TBMx    ;Output.cum.DAM.TBMx(end) +  Opt{i}.Output.cum.DAM.TBMx       ];
        Output.cum.DAM.TBMy      = [Output.cum.DAM.TBMy    ;Output.cum.DAM.TBMy(end) +  Opt{i}.Output.cum.DAM.TBMy       ];
        Output.cum.DAM.TBMz      = [Output.cum.DAM.TBMz    ;Output.cum.DAM.TBMz(end) +  Opt{i}.Output.cum.DAM.TBMz       ];
        Output.cum.DAM.BRMx      = [Output.cum.DAM.BRMx    ;Output.cum.DAM.BRMx(end) +  Opt{i}.Output.cum.DAM.BRMx       ];
        Output.cum.DAM.BRMy      = [Output.cum.DAM.BRMy    ;Output.cum.DAM.BRMy(end) +  Opt{i}.Output.cum.DAM.BRMy       ];
        Output.cum.DAM.BRMz      = [Output.cum.DAM.BRMz    ;Output.cum.DAM.BRMz(end) +  Opt{i}.Output.cum.DAM.BRMz       ];
        Output.cum.DAM.BROop     = [Output.cum.DAM.BROop   ;Output.cum.DAM.BROop(end) +  Opt{i}.Output.cum.DAM.BROop      ];
        Output.cum.DAM.BRIp      = [Output.cum.DAM.BRIp    ;Output.cum.DAM.BRIp(end) +  Opt{i}.Output.cum.DAM.BRIp       ];
        Output.cum.DAM.TTMx      = [Output.cum.DAM.TTMx    ;Output.cum.DAM.TTMx(end) +  Opt{i}.Output.cum.DAM.TTMx       ];
        Output.cum.DAM.TTMy      = [Output.cum.DAM.TTMy    ;Output.cum.DAM.TTMy(end) +  Opt{i}.Output.cum.DAM.TTMy       ];
        Output.cum.DAM.TTMz      = [Output.cum.DAM.TTMz    ;Output.cum.DAM.TTMz(end) +  Opt{i}.Output.cum.DAM.TTMz       ];
        Output.cum.DAM.LSSMy     = [Output.cum.DAM.LSSMy   ;Output.cum.DAM.LSSMy(end) +  Opt{i}.Output.cum.DAM.LSSMy      ];
        Output.cum.DAM.LSSMz     = [Output.cum.DAM.LSSMz   ;Output.cum.DAM.LSSMy(end) +  Opt{i}.Output.cum.DAM.LSSMz      ];
        Output.cum.DAM.LSSTq     = [Output.cum.DAM.LSSTq   ;Output.cum.DAM.LSSTq(end) +  Opt{i}.Output.cum.DAM.LSSTq      ];


    end


end

Output.metrics.CapacityFactor = Output.cum.Energy(end)/(10^4*length(Output.Time{1,2}));
Output.metrics.IBCactivation.hours = sum(Output.inst.IBC(:)==1);
Output.metrics.IBCactivation.Perc = 100*sum(Output.inst.IBC(:)==1)/length(Output.Time{1,2});
Output.metrics.ShutDown.hours = sum(Output.inst.Prat(:)==0);
Output.metrics.ShutDown.Perc = sum(Output.inst.Prat(:)==0)/length(Output.Time{1,2});

Output.metrics.CompToBaseCum = CompareToBaselineCum(Output,Baseline,length(Output.Time{1,2}));
Output.options.files.BaselineResult =  [pathbase BaseFile{1}]  ;

%% Functions
function CompToBase = CompareToBaselineCum(Output,Baseline,time_cnt)
CompToBase.Energy  = 100*(Output.cum.Energy(end)/Baseline.Output.cum.Energy(end)-1);
CompToBase.Revenue = 100*(Output.cum.Revenue(end)/Baseline.Output.cum.Revenue(end)-1);
CompToBase.TBMx = 100*(Output.cum.DAM.TBMx(end)/Baseline.Output.cum.DAM.TBMx(end)-1);
CompToBase.TBMy = 100*(Output.cum.DAM.TBMy(end)/Baseline.Output.cum.DAM.TBMy(end)-1);
CompToBase.TBMz = 100*(Output.cum.DAM.TBMz(end)/Baseline.Output.cum.DAM.TBMz(end)-1);
CompToBase.BRMx = 100*(Output.cum.DAM.BRMx(end)/Baseline.Output.cum.DAM.BRMx(end)-1);
CompToBase.BRMy = 100*(Output.cum.DAM.BRMy(end)/Baseline.Output.cum.DAM.BRMy(end)-1);
CompToBase.BRMz = 100*(Output.cum.DAM.BRMz(end)/Baseline.Output.cum.DAM.BRMz(end)-1);
CompToBase.BROop = 100*(Output.cum.DAM.BROop(end)/Baseline.Output.cum.DAM.BROop(end)-1);
CompToBase.BRIp = 100*(Output.cum.DAM.BRIp(end)/Baseline.Output.cum.DAM.BRIp(end)-1);
CompToBase.TTMx = 100*(Output.cum.DAM.TTMx(end)/Baseline.Output.cum.DAM.TTMx(end)-1);
CompToBase.TTMy = 100*(Output.cum.DAM.TTMy(end)/Baseline.Output.cum.DAM.TTMy(end)-1);
CompToBase.TTMz = 100*(Output.cum.DAM.TTMz(end)/Baseline.Output.cum.DAM.TTMz(end)-1);
CompToBase.LSSMy = 100*(Output.cum.DAM.LSSMy(end)/Baseline.Output.cum.DAM.LSSMy(end)-1);
CompToBase.LSSMz = 100*(Output.cum.DAM.LSSMz(end)/Baseline.Output.cum.DAM.LSSMz(end)-1);
CompToBase.LSSTq = 100*(Output.cum.DAM.LSSTq(end)/Baseline.Output.cum.DAM.LSSTq(end)-1);
CompToBase.BlPitchTrav = 100*(Output.cum.BlPitchTrav(end,1)/Baseline.Output.cum.BlPitchTrav(end,1)-1);
CompToBase.GenTqSTD    = 100*(mean(Output.inst.GenTqSTD(:,1))/mean(Baseline.Output.inst.GenTqSTD(:,1))-1);
CompToBase.GenSpeedSTD = 100*(mean(Output.inst.GenSpeedSTD(:,1))/mean(Baseline.Output.inst.GenSpeedSTD(:,1))-1);
CompToBase.BlPitchSTD  = 100*(mean(Output.inst.BlPitchSTD(:,1))/mean(Baseline.Output.inst.BlPitchSTD(:,1))-1);
CompToBase.ShutRelPerc  = 100*(Output.metrics.ShutDown.Perc/Baseline.Output.metrics.ShutDown.Perc-1);
CompToBase.ShutHoursDiff = Output.metrics.ShutDown.hours-Baseline.Output.metrics.ShutDown.hours;
CompToBase.ShutAbsHoursPerc = 100*CompToBase.ShutHoursDiff/time_cnt;
end