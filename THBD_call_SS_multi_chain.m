%% Data Loading and Setup
% P and D (both column vectors), and time vector (negative to positive)
% should be loaded here.

% modname = 'gauss2';
% data = load(['/Users/evanzhang/Library/CloudStorage/' ...
%     'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/SynMats/' ...
%     modname '.mat']);
% 
% P = data.P; % P(arent) is the reference phase: Hilbert of S
% D = data.D_model; % D(aughter) is the observation: SS
% time = data.time; % P, D, and time should have the same length

% Data (Match Filtered)
modname = 'NoMelt_matchfiltered_attenuated_1_FOM_min20';
data = load(['/Users/evanzhang/Library/CloudStorage/' ...
    'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/DataMats/' ...
    'SS_matchfiltered_NoMelt.mat']);
% P = data.S_seg_corr_avg(1:end-1);
P = data.P_attenuated;
D = data.SS_seg_corr_avg(1:end-1);
time = data.seg_time(1:end-1);

%% Input arguments
% total chain number for multi-chain runs
totalChains = 96;

% mandatory arguments
totalSteps = 10000; % total steps for MCMC
saveSteps  = 20; % save results per this many steps
output_time = 45; % half, in seconds, should not exceed input time

% check output time
if output_time > max(time)
    error('User specified output time is longer than input time.')
end

% optional arguments
% corner frequencies for optional zero-phase bandpass filtering
isBP = 0;
lf = 0.01;
hf = 1.0;

% for synthetics, filename for true solution to plot (overlay) on the model
% output; should be a *.mat file containing a field "G_model"
isPlotTrue = 0;
trueG_filename = ['/Users/evanzhang/Library/CloudStorage' ...
    '/GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/SynMats/' ...
    modname '.mat'];

% put arguments into structures
ManArgs.totalSteps  = totalSteps;
ManArgs.saveSteps   = saveSteps;
ManArgs.output_time = output_time;

OptArgs.isBP = isBP;
OptArgs.lf   = lf;
OptArgs.hf   = hf;
OptArgs.isPlotTrue = isPlotTrue;
OptArgs.trueG_filename = trueG_filename;

%% Main Call
% setup save directory
resdir = ['/Users/evanzhang/Library/CloudStorage/' ...
    'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/Results/'];
savedir = [resdir modname '/'];
if ~exist(savedir, 'dir')
    mkdir(savedir);
end

tic
parfor ichain = 1:totalChains
    
    % skip if result file exists
    savename = [savedir modname '_chain' num2str(ichain, '%03d') '.mat'];
    if exist(savename, 'file')
        continue;
    end

    % main call
    [LocSave, AmpSave, WidSave, SigSave, NC1Save, NC2Save, LikSave, prior] = ...
        THBD_SS(time, P, D, ManArgs, OptArgs); %%%%%% _Transdimensional
    
    % save results
    save_chain_result(savename, LocSave, AmpSave, WidSave, SigSave, ...
        NC1Save, NC2Save, LikSave, prior);

end
toc

%% Plot

% load all saved file names
allMat = dir([savedir '*mat']);

% load prior from the first save
savename = [savedir modname '_chain' num2str(1, '%03d') '.mat'];
tmp = load([savedir allMat(1).name]);
prior = tmp.prior;

output_time_vector = linspace(-output_time, output_time, ...
    round(2 * output_time / (time(2) - time(1))) + 1);
model_ensemble = zeros(prior.tlen * 2 - 1, length(allMat));
time_ensemble = repmat(output_time_vector(:), 1, length(allMat));

all_location = zeros(1, length(allMat));

for i = 1:length(allMat)
    savename = [savedir allMat(i).name];
    tmp = load(savename);
    model.location = tmp.LocSave(end, :);
    model.amplitude = tmp.AmpSave(end, :);
    model.width = tmp.WidSave(end, :);
    G_model = create_G_from_model_force_oceanic_Moho(model, prior); %%%%%% _transdimensional
    model_ensemble(:, i) = G_model';

    all_location(i) = model.location * prior.dt;
end

figure(1);
clf;
[N, C] = hist3([time_ensemble(:) model_ensemble(:)],[length(output_time_vector) 500]);
pcolor(C{1},C{2},N.');
shading flat; 
clim([0 size(time_ensemble, 2)/10]);
colormap(hot); 
xlabel('Time (s)');
colorbar;
title(['Model Ensemble (Average of ' num2str(length(allMat)) ' Chains)']);