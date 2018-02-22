function motion_correction(input_filenames, saving_filename)
% This function  preformes concatenatination, motion correction and cropping
% for the neuronal videos trough mosaic. It is suitible for Bambi movies,
% but trough changing filenames can work with other movies found in the
% same folder.

%% define motion correction parameters:
motionCorrectionType = 'Translation'; % 'Translation' 'Rigid','Similarity','Skew'
parallel = true; % To run computation on multiple processors
invIm = true; % Invert image black-white (imcomplement)
applyMeanFilter = true;
meanSubtraction = true;

%% Load movies
mosaic_movies = cell(1, length(input_filenames));
for i=1:length(input_filenames)
    mosaic_movies{i} = mosaic.loadMovieTiff(input_filenames{i});
end

%% Concatenate movies
movies_list = mosaic.List('mosaic.Movie', mosaic_movies);
mosaic_movie = mosaic.concatenateMovies(movies_list, 'gapTime', 0);

%% Choose roi for motion correction
msgbox('close figure for choosinf ROI');
tissueRoi = mosaic.EllipseRoi(mosaic.Point(200, 110), mosaic.Point(50, 70));
tissueRoi.edit(mosaic_movie);

%% Perform motion correction
[tifMovie_Registered ,~] = mosaic.motionCorrectMovie(...
    mosaic_movie, 'parallelProcess',parallel, ...
    'motionType',motionCorrectionType,'invertImage',invIm,'roi',tissueRoi, ...
    'subtractSpatialMean', meanSubtraction,'applySpatialMean', applyMeanFilter);

%% Crop movie
tifMovie_Registered.view()
frame1 = mosaic.extractFrame(tifMovie_Registered, 'frame', 1);
msgbox('Dubble click the ROI for chossing area for cropping')
frame1.view()
h = imrect;
position = round(wait(h));

cropped_movie = mosaic.cropMovie(tifMovie_Registered, position(1),position(2),position(3)+position(1),position(4)+position(2));
cropped_movie.view()

%% Save results
mosaic.saveMovieTiff(cropped_movie,saving_filename);
end
