% This code is to test Metal artifact simulation.
%
% Mitsuki Sakamoto <sakamoto.mitsuki.si2@is.naist.jp>
% Nara Institute of Science and Technology
% 2019-06-25
%%
addpath('.\src');
addpath('.\utils');

%% Load images
load('./sample_data/sample_2.mat'); % load "sample" valuable
image = sample.image;
metal = sample.metal;
pixel_size = sample.pixel_size; % [cm]
disp('Load images');

%%

figure(1);clf;
imagesc(image);
colorbar;
colormap(gray);

figure(2);clf;
imagesc(metal);
colorbar;
colormap(gray);

%% Set config
config = set_config_for_artifact_simulation(pixel_size);
disp('Set config');

%% Preprocess
image(image<-500) = -1000; % erase the boundary
image = hu2mu(double(image), config.mu_water, config.mu_air);
disp('Preprocess');

%%

figure(3);clf;
imagesc(image);
colorbar;
colormap(gray);

%% Phantom calibration

phantom = create_phantom(512, 512, 200, config.mu_water);
config.correction_coeff = water_correction(phantom, config);
disp('Phantom calibration');

%%

figure(4);clf;
imagesc(phantom);
colorbar;
colormap(gray);

%% Metal Artifact Simulation
sim = metal_artifact_simulation(image, metal, config);
disp('Metal Artifact Simulation');

%%

figure(4);clf;
imagesc(sim);
colorbar;
colormap(gray);


%% Convert results from mu to HU
sim_hu = mu2hu(sim, config.mu_water, config.mu_air);
disp('Convert results from mu to HU');

%%

figure(5);clf;
imagesc(sim_hu);
colorbar;
colormap(gray);
caxis([0, 0.1*max(sim_hu(:))]);

%% Save results
save_dir = './outputs';
if ~exist(save_dir, 'dir'); mkdir(save_dir); end
imwrite(set_window(mu2hu(image, config.mu_water, config.mu_air), -150, 350),...
        fullfile(save_dir, 'input.png'));
imwrite(set_window(sim_hu, -150, 350),...
        fullfile(save_dir, 'output.png'));

if ~verLessThan('matlab', '9.1') % older than 2016b
  save_config_as_json(fullfile(save_dir, 'simulation_config.json'), config);
else
  save(fullfile(save_dir, 'simulation_config.mat'), 'config')
end
disp('Save results');

%%
x_metal = metal;
% parse arguments
data = config.data;
energy_composition = config.energy_composition;
E0 = config.E0;
mu_air = config.mu_air;
metal_name = config.metal_name;
metal_density = config.metal_density;
T1 = config.T1;
T2 = config.T2;
noise_scale = config.noise_scale;
filter_name = config.filter_name;
freqscale = config.freqscale;
correction_coeff = config.correction_coeff;

SOD = config.SOD;
angle_size = config.angle_size;
angle_num = config.angle_num;
pixel_size = config.pixel_size;
output_size = config.output_size;

m0_water = data{E0, 'Water'};
m0_bone = data{E0, 'Bone'};
m0_metal = data{E0, metal_name};
mu_water0 = m0_water * 1.0;
mu_metal0 = m0_metal * metal_density;

% Threshold-based weighting
T1 = hu2mu(T1, mu_water0, mu_air);
T2 = hu2mu(T2, mu_water0, mu_air);
[x_water, x_bone] = threshold_based_weighting(image, T1, T2);
x_water(x_metal>0) = 0;
x_bone(x_metal>0) = 0;
x_metal = double(x_metal) * mu_metal0;


figure(1);clf;
imagesc([x_water, x_bone, x_metal]);
colorbar;
colormap(gray);

%%


% Forward Projection
d_water = fanbeam(x_water, ...
    SOD,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing', angle_size, ...
    'FanRotationIncrement',360/angle_num);
d_bone = fanbeam(x_bone,...
    SOD,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing', angle_size, ...
    'FanRotationIncrement',360/angle_num);
d_metal = fanbeam(x_metal,...
    SOD,...
    'FanSensorGeometry','arc',...
    'FanSensorSpacing', angle_size, ...
    'FanRotationIncrement',360/angle_num);
d_water = d_water * pixel_size;
d_bone = d_bone * pixel_size;
d_metal = d_metal * pixel_size;

% Energy Composition
total_intensity = 0;
v = zeros(size(d_water, 1), size(d_water, 2), numel(energy_composition));

for ii = 1:numel(energy_composition)
    energy = energy_composition(ii);
    
    m_water = data{energy, 'Water'};
    m_bone = data{energy, 'Bone'};
    m_metal = data{energy, metal_name};
    intensity = data{energy, 'Intensity'};
    d_water_tmp = d_water*(m_water/m0_water);
    d_bone_tmp = d_bone*(m_bone/m0_bone);
    d_metal_tmp = d_metal*(m_metal/m0_metal);
    DRR = d_water_tmp + d_bone_tmp + d_metal_tmp;
    
    y = intensity * (exp(-DRR));
    v(:, :, ii) = y;
    total_intensity = total_intensity + intensity;
end
poly_y = sum(v, 3);

figure(2);clf;
imagesc(poly_y);
colorbar;
colormap(gray);

p = -log(poly_y*0.0000001); 
sim = ifanbeam(p,...
               SOD,...
               'FanSensorGeometry','arc',...
               'FanSensorSpacing',angle_size,...
               'OutputSize',output_size,... 
               'FanRotationIncrement',360/angle_num,...
               'Filter', filter_name,...
               'FrequencyScaling', freqscale);
sim(sim<0) = 0;
sim = sim/pixel_size;

figure(3);clf;
imagesc(sim);
colorbar;
colormap(gray);

caxis([0, 0.25*max(sim(:))]);







































