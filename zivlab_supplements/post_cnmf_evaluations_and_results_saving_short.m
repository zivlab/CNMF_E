
prompt = {'microns per pixel'};

dlg_title = 'please insert the following details';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

%mouse_name=['C' answer{1} 'M' answer{2}]; %         try again- C67M3- 2.55  ; 3 10- C107m3- 2.6  % ?67?3- ???????

%session=answer{3};
microns_per_pixel=str2double(answer{1});

%% creating a stamp for the relevant parameters set
gSig=neuron.options.gSig;
gSiz=neuron.options.gSiz;
min_corr=neuron.options.min_corr;
min_pnr=neuron.options.min_pnr;

parameters_stamp=[num2str(gSig) ' ' num2str(gSiz) ' ' num2str(min_corr) ' ' num2str(min_pnr)];

%% neurons display option
neurons_display_choice = menu('would you like to save the filters and traces visualization? ','Yes','No');

%% defining the folder in which the results wil be saved
h = msgbox('please select the folder in which the results will be saved');
pause(2)
close all hidden;
start_path='Z:\Short term data storage\Lab members';
dialog_title='results folder';
results_path = uigetdir(start_path,dialog_title);


mkdir(results_path,  parameters_stamp);

current_results_path=[results_path,  parameters_stamp];


%% display contours of the neurons
Cn=neuron.Cn;

neuron.show_contours(0.6);
colormap gray;
title(['contours of estimated neurons (' num2str(size(neuron.A,2)) ' neurons)']);

saveas(gcf,[current_results_path 'all_contours.fig']);
saveas(gcf,[current_results_path 'all_contours.png']);

%% creating correlation and PNR image

[cn, pnr] = neuron.correlation_pnr_parallel([1, 5000]);

tmp_d = max(1,round(gSiz/4));
v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
ind = (v_max==cn.*pnr);

correlation_and_pnr_figure=figure('papersize', [12, 3]);
init_fig;
subplot(131);
imagesc(cn, [0, 1]);
title('local corr. image');
axis equal off tight;
subplot(132);
pnr_vmax = max(pnr(:))*0.8;
imagesc(pnr, [3, pnr_vmax]);
axis equal off tight;
title('PNR image');
subplot(133);
imagesc(cn.*pnr, [3, pnr_vmax]);
hold on;

tmp_ind = ind & (cn>=min_corr) & (pnr>=min_pnr);
[r, c] = find(tmp_ind);
ax_seeds = plot(c, r, '.m', 'markersize', 10);
axis equal off tight;
title('candidate seed pixels');
ylabel('PNR');
xlabel('Cn');tmp_d = max(1,round(gSiz/4));
v_max = ordfilt2(cn.*pnr, tmp_d^2, true(tmp_d));
ind = (v_max==cn.*pnr);

%% reshaping spatial filters
d1=neuron.options.d1;
d2=neuron.options.d2;

number_of_cells=size(neuron.A,2);

allFiltersMat=nan(number_of_cells,d1,d2);

for n=1:number_of_cells
    allFiltersMat(n,:,:)=reshape(neuron.A(:,n), d1, d2);
    
end

%%

%% Shape sorting from preprocessing GUI

edges_flag=false;
max_radius=30;
min_radius=10;

max_area=pi*max_radius^2;
min_area=pi*min_radius^2;
options.maxCellSize=round(max_area/microns_per_pixel^2); % Maximum area of cell in region props
options.minCellSize=round(min_area/microns_per_pixel^2); %  Minumum area of cell in region props

options.circularityThreshold=0.8;
options.secRoiThrshld = 0.1; % maximum number of pixels that we can drop from the filter  and still keep it as a cell (if region props finds 2 rois)

sizeOfImage = size(allFiltersMat);
sizeOfImage=sizeOfImage(2:3);
edgePixels = GetLinearIndicesOfImageEdge(sizeOfImage);

filtersFeatures = cell(number_of_cells,1);
idxOfGoodICs = [];

circularity=zeros(1,number_of_cells);
h = waitbar(0,'Sorting cells according to typical cellular morphology');
for n=1:1:number_of_cells
    waitbar((n)/number_of_cells,h,['Shape sorting - inspecting cell number ' num2str(n) '/' num2str(number_of_cells)])
    [filtersFeatures{n},probICFlag,circularity(n),distance(n)] = Check_Features(squeeze(allFiltersMat(n,:,:)),edgePixels,options,edges_flag);
    %     probICFlag=0;
    if ~probICFlag
        idxOfGoodICs = cat(1,idxOfGoodICs,n);
    end
end
close(h);

all_area=cellfun (@(x) x.Area,filtersFeatures,'UniformOutput',false);
all_area = cell2mat(all_area);
all_radius=sqrt(all_area./pi).*microns_per_pixel;
center_x=round(size(allFiltersMat,3)/2);
center_y=round(size(allFiltersMat,2)/2);
max_distance=max(distance);
allowed_max_area=[options.maxCellSize options.maxCellSize+options.maxCellSize];
allowed_min_area=[options.minCellSize options.minCellSize+options.minCellSize];
allowed_circularity=[options.circularityThreshold options.circularityThreshold-options.circularityThreshold/3];
max_radius=sqrt(allowed_max_area./pi).*microns_per_pixel;
min_radius=sqrt(allowed_min_area./pi).*microns_per_pixel;
distance=distance.*microns_per_pixel;


%% by pp gui size & circularity sorting - insert as output for the function

shape_sorting_excluded_neurons=~ismember(1:size(neuron.A,2),idxOfGoodICs);
shape_sorting_excluded_neurons_idx=find(shape_sorting_excluded_neurons);


num_excluded_neurons=length(shape_sorting_excluded_neurons_idx);
%Cn = imresize(Cn, [d1, d2]);
neuron.show_contours(0.6);
colormap gray;

for i=1:num_excluded_neurons
    
    chosen_neuron=neuron.A(:,shape_sorting_excluded_neurons_idx(i));
    
    %     Cn = imresize(Cn, [d1, d2]);
    %     neuron.show_contours(0.6);
    %     colormap gray;
    
    
    
    % add also a visualization of the chosen neuron
    hold on
    
    % figure;
    thr=0.9;
    A_temp = full(reshape(chosen_neuron,d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    
    contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'linecolor',[0 0.5 0.5],'linewidth', 3);
    
    
    
    % Cn = imresize(Cn, [d1, d2]);
    % neuron.show_contours(0.6);
    % colormap gray;
    
    % add also a visualization of the chosen neuron
    hold on
    
    % figure;
    
    
    
    
end
title(['contours of shape sorted neurons (' num2str(length(idxOfGoodICs)) ' neurons)']);

shape_sorting_figure=figure;
subplot(1,2,1)
scatter(distance/microns_per_pixel,all_radius,100,'.')
hold on
plot([0 max_distance],min_radius,'color','r','linewidth',2)
hold on
plot([0 max_distance],max_radius,'color','r','linewidth',2)

hold on
scatter(distance(shape_sorting_excluded_neurons_idx)/microns_per_pixel,all_radius(shape_sorting_excluded_neurons_idx),100,'r.')


ylim([0 50])
xlabel('Distance from center of FOV [\mum]','FontWeight','bold')
ylabel('Cell radius [\mum]','FontWeight','bold')
subplot(1,2,2)
scatter(distance/microns_per_pixel,circularity,100,'.')
hold on
plot([0 max_distance],allowed_circularity,'color','r','linewidth',2)

hold on

scatter(distance(shape_sorting_excluded_neurons_idx)/microns_per_pixel,circularity(shape_sorting_excluded_neurons_idx),100,'r.')






ylim([0 2])
xlabel('Distance from center of FOV [\mum]','FontWeight','bold')
ylabel('Circularity','FontWeight','bold')


%% calculating neurons' centers

centers=zeros(n,2);

for n=1:number_of_cells
    
    centers(n,:) = filtersFeatures{n}.Centroid;
    
end


%% sorting results matrices by shapes sorting

finalFiltersMat=allFiltersMat(idxOfGoodICs,:,:);
finalTracesMat=neuron.C(idxOfGoodICs,:);
finalSpikesMat=neuron.S(idxOfGoodICs,:);
finalNeuronsCenters=centers(idxOfGoodICs,:);
%%  looking at number of joint pixels between neurons-
% for every cell (column)- the joint pixels fraction is indicated by the entry along the rows. every row- is for another cell

neurons_locations=neuron.A;
neurons_locations(neurons_locations~=0)=1;
sum_of_joint_pixels = (neurons_locations)' *(neurons_locations);
sum_of_pixels_per_neuron=sum(neurons_locations);
rep_sum_of_pixels_per_neuron=repmat(sum_of_pixels_per_neuron,size(sum_of_joint_pixels,1),1);
joint_pixels_measure=sum_of_joint_pixels./rep_sum_of_pixels_per_neuron; % for each neuron- the fraction of his pixels are in cluded in other contours appears aling its column
joint_pixels_measure_with_others=joint_pixels_measure;
joint_pixels_measure_with_others(~~eye(size(joint_pixels_measure)))=NaN;


joint_pixels_figure=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
histogram(joint_pixels_measure_with_others(joint_pixels_measure_with_others>0),0:0.01:1,'normalization','probability');
ylabel('number of cells');
xlabel('overlap');
ylim([0 0.5]);

title('joint pixels distribution (given overlap>0)');

subplot(1,2,2)
histogram(joint_pixels_measure_with_others(joint_pixels_measure_with_others>0.5),0.7:0.01:1,'normalization','probability');
ylabel('number of cells');
xlabel('overlap');
ylim([0 0.4]);

title('joint pixels distribution(given overlap>0.7)');
suptitle('joint pixels estimation');
%% SNR and noise

SNR_post_updates = var(neuron.C, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);

SNR_post_updates_c_raw = var(neuron.C_raw, 0, 2)./var(neuron.C_raw-neuron.C, 0, 2);

noise_estimation=var(neuron.C_raw-neuron.C, 0, 2);

SNR_figure=figure;
histogram(SNR_post_updates,0:0.25:20,'normalization','probability');
ylim([0 0.3]);
xlabel('SNR');
ylabel('probability');


%% results saving
save([current_results_path 'finalFiltersMat.mat'],'finalFiltersMat','-v7.3');
save([current_results_path 'finalTracesMat.mat'],'finalTracesMat','-v7.3');
save([current_results_path 'finalSpikesMat.mat'],'finalSpikesMat','-v7.3');
save([current_results_path 'finalNeuronsCenters.mat'],'finalNeuronsCenters','-v7.3');



results.A=neuron.A;
results.C=neuron.C;
results.S=neuron.S;
results.C_raw=neuron.C_raw;
results.allFiltersMat=allFiltersMat;
results.allCenters=centers;

results.microns_per_pixel=microns_per_pixel;
results.idxOfGoodICs=idxOfGoodICs;

results.Cn=neuron.Cn;
results.P=neuron.P;
results.options=neuron.options;

results.video_file_name=nam;


parameters.with_dendrites=with_dendrites;

results.merge_thr_spatial=merge_thr_spatial;

results.min_pixel=min_pixel;
results.min_pnr=min_pnr;
results.min_corr=min_corr;
%results.Nspatial=Nspatial;
%results.maxIter=maxIter;
%results.additional_search_flag=additional_search_flag;

results.shape_sorting_options=options;

% estimates per neuron:
results.joint_pixels_measure_with_others=joint_pixels_measure_with_others; % Noa's measure
results.noise_estimation=noise_estimation;
results.SNR_post_updates=SNR_post_updates;


save([current_results_path 'results.mat'],'results','-v7.3');

%% with hadas- choosing spesific figures

saveas(correlation_and_pnr_figure, [current_results_path 'correlation_and_pnr.fig']);
saveas(correlation_and_pnr_figure, [current_results_path 'correlation_and_pnr.png']);

saveas(shape_sorting_figure, [current_results_path 'shape_sorting.fig']);
saveas(shape_sorting_figure, [current_results_path 'shape_sorting.png']);

saveas(joint_pixels_figure, [current_results_path 'joint_pixels.fig']);
saveas(joint_pixels_figure, [current_results_path 'joint_pixels.png']);

saveas(SNR_figure, [current_results_path 'SNR.fig']);
saveas(SNR_figure, [current_results_path 'SNR.png']);



%% display neurons

if neurons_display_choice
    mkdir(current_results_path, ['neurons- session ' session]);
    
    current_neurons_path=[current_results_path '\neurons'];
    
    dir_neurons = sprintf('%s%s%s_neurons%s', current_neurons_path);
    neuron.save_neurons(current_neurons_path);
end


%% create a video for displaying the
% amp_ac = 140;
% range_ac = 5+[0, amp_ac];
% multi_factor = 10;
% range_Y = 1300+[0, amp_ac*multi_factor];
% 
% avi_filename = neuron.show_demixed_video(save_demixed, kt, [], amp_ac, range_ac, range_Y, multi_factor);

%%

msgbox('finished saving results');