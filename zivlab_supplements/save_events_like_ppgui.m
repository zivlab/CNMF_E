% This script runs events detection on traces mat outside the PPGUI
traces_filename = 'Z:\experiments\projects\bambi\linear_track_2\analysis\2_groups\c51m4\day_1\CNMF-E\smoothed_neuronal_source_extraction\frames_1_9000\LOGS_20-Feb_11_18_43\5 8 0.8 8\finalTracesMat.mat';
eventParams.eventTh = 5;
eventParams.tau = 0.2;
eventParams.fs = 10;
eventParams.rateTh = 0.01;

load(traces_filename);
traces = finalTracesMat';

[eventsMatrix,onsetWidthMatrix]=event_detection(traces,eventParams);
[idxOfGoodICs] = trace_sorting(traces,eventParams.rateTh,eventParams.tau,eventParams.fs);
%%
allEventsMat = eventsMatrix(:, idxOfGoodICs)';
allTracesMat = finalTracesMat;
save('finalEventsMat.mat', 'allEventsMat')
save('finalTracesMat.mat', 'allTracesMat')