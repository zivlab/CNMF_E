% Prepering movies for CNMFE head fix experiment
data_path = 'E:\experiments\head_fix_bambi\data\day7\C81M6';
list_dirs = dir(data_path);
list_dirs = list_dirs(3:end);
filenames = {};

%%
j=1
for i=1:length(list_dirs)
    if list_dirs(i).isdir
        name_length = length(list_dirs(i).name);
        date{j} = list_dirs(i).name(name_length - 10:end);
        dir_name{j} = list_dirs(i).name;
        j=j+1;
    end
end

[~, ind] = sort(date);
dir_name = dir_name(ind);

for i=1:length(dir_name)
%     filenames{2*i-1} = [data_path, '\', dir_name{i}, '\f0_frames.tif'];
    filenames{i} = [data_path, '\', dir_name{i}, '\dfof_frames.tif'];
   
end

%%
mosaic.initialize;
saving_filename = [data_path, '\concatenated_neuronal.tif'];
motion_correction(filenames, saving_filename, ...
                           1, 0)
save([data_path, '\recording_list.mat'], 'filenames');
                       
mosaic.terminate;               