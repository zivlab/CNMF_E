clear all; close all; clc;
EVENT_RATE_TRESH = 0.01;
SNR_THRESH = 0.5;

analysis_path = 'D:\experiments_processing\head_fix_treadmill_bambi\c60m1';
files = dir(analysis_path);
for i=1:length(files)
    if strfind(files(i).name, 'day_0')
        %% Sort
        results_path = [analysis_path, '\', files(i).name, '\CNMFE\3 8 0.7 8\results.mat']
        load(results_path);
        spikes = results.S;
        event_rate = mean(spikes, 2);
        event_rate_ind = event_rate > EVENT_RATE_TRESH;
        
        snr = results.SNR_post_updates;
        snr_ind = snr > SNR_THRESH;
        
        number_of_neurons = length(snr);
        good_ics = zeros(number_of_neurons, 1, 'logical');
        good_ics(results.idxOfGoodICs) = true;
        nrn_sort_ind = event_rate_ind & snr_ind & good_ics;
        
        results.nrn_sort_ind = nrn_sort_ind;
        
        %% Save all the variables after sorting
        finalFiltersMat = results.allFiltersMat(nrn_sort_ind,:,:);
        finalTracesMat = results.C(nrn_sort_ind,:);
        finalSpikesMat = results.S(nrn_sort_ind,:);
        finalNeuronsCenters=results.allCenters(nrn_sort_ind,:);
        
        save_path = [analysis_path, '\', files(i).name, '\cooked2'];
        mkdir(save_path)
        save(fullfile(save_path, 'finalFiltersMat.mat'),'finalFiltersMat','-v7.3');
        save(fullfile(save_path, 'finalTracesMat.mat'),'finalTracesMat','-v7.3');
        save(fullfile(save_path, 'finalSpikesMat.mat'),'finalSpikesMat','-v7.3');
        save(fullfile(save_path, 'finalNeuronsCenters.mat'),'finalNeuronsCenters','-v7.3');
        save(fullfile(save_path, 'results.mat'),'results','-v7.3');
    end
end
    