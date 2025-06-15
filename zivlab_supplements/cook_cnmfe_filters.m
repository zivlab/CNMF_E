% cook cnmfe filters
filters_filename = "D:\experiments_processing\bambi_treadmill_04-21\day_2\c116m8\cnmfe\3 8 0.85 8\finalFiltersMat.mat";
save_path = "D:\experiments_processing\bambi_treadmill_04-21\day_2\c116m8\cnmfe\3 8 0.85 8\cookedFiltersMat.mat";
number_of_neurons = size(finalFiltersMat, 1);
cooked_filters_mat = zeros(size(finalFiltersMat));
for i=1:number_of_neurons
    curr_filter = squeeze(finalFiltersMat(i, :, :));
    threshold = 0.5 * max(max(curr_filter));
    ind = curr_filter >= threshold;
    b = zeros(size(curr_filter));
    b(ind) = curr_filter(ind);
    cooked_filters_mat(i, :, :) = b;
%     figure;
%     subplot(1, 2, 1);
%     imshow(curr_filter, []); colormap('jet');
%     subplot(1, 2, 2);
%     imshow(squeeze(cooked_filters_mat(i,:, :)), []); colormap('jet');
end

allFiltersMat = cooked_filters_mat;
save(save_path, 'allFiltersMat');