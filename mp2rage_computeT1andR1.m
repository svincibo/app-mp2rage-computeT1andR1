function [] = mp2rage_computeT1andR1()

% app-mp2rage-computeT1andR1

% Computes T1map and R1map when provided a denoised unit1.nii.gz image.

% Setup environment.
if ~isdeployed
    
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('nii_func'));
    addpath(genpath('func'));
    
end

% Read in config.json.
config = jsondecode(fileread('config.json'));

% Read in config.json files that are contain the meta-data for each image.
config_inv1 = jsondecode(fileread(config.json_inv1));
config_inv2 = jsondecode(fileread(config.json_inv2));

% Assign constants related to the mp2rage acquisition protocol.
MP2RAGE.B0 = config_inv1.MagneticFieldStrength; % in Tesla
MP2RAGE.TR = config_inv1.RepetitionTime; % MP2RAGE TR in seconds
MP2RAGE.TRFLASH = (config_inv2.InversionTime - config_inv1.InversionTime)/(config_inv1.BaseResolution*(config_inv1.PercentPhaseFOV/100)); % TR of the GRE readout, per Hu: (TI2-TI1)/(number of phase encoding steps) = (2.5-0.7)/(256*.938) ms = 7.5 ms.
MP2RAGE.TIs = [config_inv1.InversionTime config_inv2.InversionTime];% inversion times - time between middle of refocusing pulse and excitatoin of the k-space center encoding
if isempty(config.slicesperslab)
          error('Sorry -- this app is not yet able to determine the SlicesPerSlab parameter from your metadata. Please look at your MRI sequence protocol to find this value and enter it manually when running the app.');
	  MP2RAGE.NZslices=config_inv1.SlicesPerSlab.*[config_inv1.PartialFourier-0.5 0.5]; % Slices Per Slab * [PartialFourierInSlice-0.5  0.5] 
else
	MP2RAGE.NZslices=config.slicesperslab.*[config_inv1.PartialFourier-0.5 0.5]; % Slices Per Slab * [PartialFourierInSlice-0.5  0.5] ;
end
MP2RAGE.FlipDegrees = [config_inv1.FlipAngle config_inv2.FlipAngle]; % Flip angle of the two readouts in degrees

%% Calculate T1map and R1map.

% Read in the denoised UNI images, i.e., the mp2rage image.
MP2RAGE.filename = config.unit1; % standard MP2RAGE T1w image;
MP2RAGE.filenameOUT = fullfile('output', 'unit1.nii.gz');% image without background noise;

if ~isempty(config.mask)
    
    % Read in mask.
    mask = niftiRead(config.mask);
    
    % Convert mask.data class from int32 to double for compatability with MP2RAGEimg.
    mask.data = double(mask.data);

end

%% Calculate T1map and R1map.

% load the MP2RAGE data - it can be either the SIEMENS one scaled from
% 0 4095 or the standard -0.5 to 0.5
clear MP2RAGEimg;
MP2RAGEimg = load_untouch_nii(MP2RAGE.filename);

% Caluclate the T1 and R1 maps for this subject.
clear T1map R1map
[T1map, R1map] = T1estimateMP2RAGE(MP2RAGEimg, MP2RAGE, 0.96);

% Convert MP2RAGEimg.img class from int16 to double for compatability with mask.data.
MP2RAGEimg.img = double(MP2RAGEimg.img);

%% Create and write out T1 and R1 nifti files.

% Intensity values.
mp2rage_intensity = MP2RAGEimg.img;

%T1 values
nii_temp = niftiRead(config.unit1);
nii_temp.data = T1map.img;
niftiWrite(nii_temp, fullfile('output', 'T1.nii.gz'));
if ~isempty(config.json_unit1)
    copyfile(config.json_unit1, fullfile('output', 'T1.json'));
end
T1 = T1map.img;

% R1 values
nii_temp.data = R1map.img;
niftiWrite(nii_temp, fullfile('output', 'R1.nii.gz'));
if ~isempty(config.json_unit1)
    copyfile(config.json_unit1, fullfile('output', 'R1.json'));
end
R1 = R1map.img;

%% QA: Construct B1 sensitivity plot for our specific acquisition parameters.

figure(1); hold on;
plotMP2RAGEproperties(MP2RAGE)
saveas(gcf, 'qa1.png'); 
hold off;

%% QA: Plot histogram of signal intensity, T1, and R1.

figure(2); hold on;
facealpha = 0.4; edgecolor = 'none'; linewidth = 2;

% Histogram of intensity values, not including zeros
mp2rage_intensity(find(mp2rage_intensity==0)) = NaN;
temp = mp2rage_intensity(:);
mp2rage_intensity_out = temp(~isnan(temp)); clear temp;

subplot(1, 3, 1)
facecolor = [128 128 128]/255; %gray
linecolor = facecolor;
toplot = double(mp2rage_intensity_out(:));
histogram(toplot, 'Normalization', 'probability', 'EdgeColor', edgecolor, 'FaceColor', facecolor, 'FaceAlpha', facealpha);
hold on;
[N, edges] = histcounts(toplot, 'Normalization', 'probability');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N, 'Color', linecolor, 'LineWidth', linewidth);
pbaspect([1 1 1]);
ylabel('Frequency');
xlabel('Signal Intensity');
hold off;

% Histogram of T1 values, not including zeros
T1(find(T1==0)) = NaN;
temp = T1(:);
T1_out = temp(~isnan(temp)); clear temp;

subplot(1, 3, 2)
clear toplot;
facecolor = [40 62 102]/255; %dark blue
linecolor = facecolor;
toplot = T1_out(:);
histogram(toplot, 'Normalization', 'probability', 'EdgeColor', edgecolor, 'FaceColor', facecolor, 'FaceAlpha', facealpha);
hold on;
[N, edges] = histcounts(toplot, 'Normalization', 'probability');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N, 'Color', linecolor, 'LineWidth', linewidth);
pbaspect([1 1 1]);
ylabel('Frequency');
xlabel('T1, seconds');
hold off;

% Histogram of R1 values, not including zeros
R1(find(R1==0)) = NaN;
temp = R1(:);
R1_out = temp(~isnan(temp)); clear temp;

subplot(1, 3, 3)
clear toplot;
facecolor = [104 149 197]/255; %light blue
linecolor = facecolor;
toplot = R1_out(:);
histogram(toplot, 'Normalization', 'probability', 'EdgeColor', edgecolor, 'FaceColor', facecolor, 'FaceAlpha', facealpha);
hold on;
[N, edges] = histcounts(toplot, 'Normalization', 'probability');
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges, N, 'Color', linecolor, 'LineWidth', linewidth);
pbaspect([1 1 1]);
ylabel('Frequency');
xlabel('R1, (seconds^-^1)');

saveas(gcf, 'qa2.png'); 
hold off;
