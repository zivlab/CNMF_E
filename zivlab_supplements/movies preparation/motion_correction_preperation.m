PATH_DIR = 'E:\linor\triangular pilot\C95M3\DAY3';
saving_filename = 'E:\linor\triangular pilot\C95M3\DAY3\concatenated.tif';
filenames = dir(PATH_DIR);

trials_filenames = {};
j=1;
for i=1:length(filenames)
    f = filenames(i).name;
    if isempty(strfind(f,'recording'))
        continue;
    end
    full_path = [PATH_DIR, '\', f];
    trials_filenames{j} = full_path;  
    j = j + 1;
end
mosaic.initialize()
motion_correction(trials_filenames, saving_filename, 1, 1)
mosaic.terminate()