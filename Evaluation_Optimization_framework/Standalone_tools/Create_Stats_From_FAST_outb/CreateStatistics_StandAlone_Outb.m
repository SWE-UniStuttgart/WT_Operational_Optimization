function Stats = CreateStatistics_StandAlone_Outb(TimeToCut,ResultFileName)

% Stand alone function if something fails and need to collect the files and
% post process the stats again
%
% Dependency: sig2ext rainflow ReadFASTbinary
%
% -ResultFileName: full path to the file including file name not the stats part. 
% -Time to cut time in seconds to cut from the simulation
%
% Example: Stats =CreateStatistics_StandAlone_Outb(60,'ExampleFiles\constTSR_IPC_P8d5_WSP18_TI02_SD1.SFunc.outb'
%
% Vasilis Pettas, Stuttgart Wind Energy (SWE), University of Stuttgart

%% Create inputs

[DataAll_init,Chanels,~,~,~] = ReadFASTbinary(ResultFileName);
TimeVec = round([DataAll_init(1,1):DataAll_init(end,1)/(length( DataAll_init(:,1))-1):DataAll_init(end,1)]',3);
timeinterval = find(TimeVec(:,1)==TimeToCut):size(TimeVec,1);
TimeVec = TimeVec(timeinterval,:);
DataAll = DataAll_init(timeinterval,:);


%% Do calculations for all the channels (first channel is time!)

Statistics.mean   =  mean(DataAll);
Statistics.median =  median(DataAll);
Statistics.min  = min(DataAll);
Statistics.max  = max(DataAll);
Statistics.std  = std(DataAll);
Statistics.MaxAbs = max(abs(DataAll));
Statistics.MinAbs = min(abs(DataAll));
Statistics.Sy    = sum(DataAll);
Statistics.SyAbs = sum(abs(DataAll));


dt =  TimeVec(2,1)-TimeVec(1,1);
N_REF = TimeVec(end,1)-TimeVec(1,1); % 1Hz DEL
vWindow    = hamming(floor((length(DataAll(:,1)))/36/2)*2);  % orig 8
for i=1:size(DataAll,2)
    [ext, exttime] = sig2ext(double(1e3*DataAll(:,i)),dt);
    rf = rainflow(ext,exttime);
    S{i} = rf(1,:); % Cycles amplitude,
    N{i} = rf(3,:); % Number of cycles
    DEL1Hz_4_in(i)  = sum((S{i} .^4).*(N{i}/N_REF)).^(1/4);
    DEL1Hz_10_in(i) = sum((S{i} .^10).*(N{i}/N_REF)).^(1/10); %#ok<*AGROW>
    [PSD_vec{i},f_vec{i}] = pwelch(detrend(DataAll(:,i),'constant'),vWindow,[],[],1/dt,'onesided'); 
    clear ext exttime rf 
end
Row_names = Chanels;
Mean = Statistics.mean';
Min = Statistics.min';
Max = Statistics.max';
Std = Statistics.std';
MaxAbs = Statistics.MaxAbs';
MinAbs = Statistics.MinAbs';
Sy = Statistics.Sy';
SyAbs = Statistics.SyAbs';
DEL1Hz_4 = DEL1Hz_4_in';
DEL1Hz_10 =DEL1Hz_10_in';
S_amp = S';
N_cycl = N';
PSD = PSD_vec';
f = f_vec';
Stats = table(Mean,Max,Min,Std,MaxAbs,MinAbs,Sy,SyAbs,DEL1Hz_4,DEL1Hz_10,S_amp,N_cycl,f,PSD,'RowNames',Row_names);

%% save the variable with the correct name
save('-v7.3',[ResultFileName(1:end-11),'_results_stats.mat'],'Stats')
% could be used to send it to the correct place
% movefile(fullfile(ResultsFilesDir,[SimulationName,'_results_stats.mat']),[ResultsFilesDir,'Stats']);



