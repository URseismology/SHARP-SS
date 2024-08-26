%% Data Loading and Setup
% P and D (both column vectors), and time vector (negative to positive)
% should be loaded here.

% Synthetics
data = load(['/Users/evanzhang/Library/CloudStorage/' ...
    'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/SynMats/' ...
    'Noise/gauss2_OS/gauss2_OS_Noise01.mat']);

P = data.P; % P(arent) is the reference phase: -Hilbert of S
D = data.D_model; % D(aughter) is the observation: SS
time = data.time; % P, D, and time should have the same length

% % align at maximum and cut
% dt = time(2) - time(1);
% cut_len = 50;
% [~, maxind_P] = max(P);
% P = P(maxind_P - round(cut_len/dt):maxind_P + round(cut_len/dt));
% P = P / max(abs(P));
% [~, maxind_D] = max(D);
% D = D(maxind_D - round(cut_len/dt):maxind_D + round(cut_len/dt));
% D = D / max(abs(D));
% time = (0:1:length(P)-1) * dt - cut_len;

% % Data (Hilbert Transform)
% data = load(['/Users/evan/Library/CloudStorage/' ...
%     'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/DataMats/' ...
%     'SS_DataTest_G_ATD.mat']);
% P = -data.hilb_S';
% D = data.orig_SS';
% time = data.time;

% % Data (Ocean Stack as Reference)
% data = load(['/Users/evan/Library/CloudStorage/' ...
%     'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/DataMats/' ...
%     'HYB_avg.mat']);
% D = data.D';
% time = data.time; % P, D, and time should have the same length
% OceanStack = load(['/Users/evan/Library/CloudStorage/' ...
%     'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/DataMats/' ...
%     'OceanStack.mat']);
% P = OceanStack.SS_mean_10Hz';

% % Data (Match Filtered)
% data = load(['/Users/evanzhang/Library/CloudStorage/' ...
%     'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/DataMats/' ...
%     'SS_matchfiltered_NoMelt.mat']);
% % P = data.S_seg_corr_avg(1:end-1);
% P = data.P_attenuated;
% D = data.SS_seg_corr_avg(1:end-1);
% time = data.seg_time(1:end-1);

%% Input arguments
% mandatory arguments
totalSteps = 10000; % total steps for MCMC
saveSteps  = 5; % save results per this many steps
output_time = 45; % half, in seconds, should not exceed input time

% check output time
if output_time > max(time)
    error('User specified output time is longer than input time.')
end

% optional arguments
% corner frequencies for optional zero-phase bandpass filtering
isBP = 0;
lf = 0.005;
hf = 0.1;

% for synthetics, filename for true solution to plot (overlay) on the model
% output; should be a *.mat file containing a field "G_model"
isPlotTrue = 1;
trueG_filename = ['/Users/evanzhang/Library/CloudStorage' ...
    '/GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/SynMats' ...
    '/gauss2_OS.mat'];

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

tic
[LocSave, AmpSave, WidSave, SigSave, NC1Save, NC2Save, LikSave, prior] = ...
    THBD_SS_Transdimensional(time, P, D, ManArgs, OptArgs); %%%%%% _Transdimensional
toc


%% Plot

end_count = 100;
output_time_vector = linspace(-output_time, output_time, ...
    round(2 * output_time / (time(2) - time(1))) + 1);
model_ensemble = zeros(prior.tlen * 2 - 1, end_count);
time_ensemble = repmat(output_time_vector(:), 1, end_count);
for i = 1:end_count
    model.location = LocSave(end - i + 1, :);
    model.amplitude = AmpSave(end - i + 1, :);
    model.width = WidSave(end - i + 1, :);
    G_model = create_G_from_model_transdimensional(model, prior); %%%%%% _transdimensional
    model_ensemble(:, i) = G_model';
end

figure(2);
clf;
subplot(2,3,1:3);
[N, C] = hist3([time_ensemble(:) model_ensemble(:)],[length(output_time_vector) 100]);
pcolor(C{1},C{2},N.');
shading flat; 
clim([0 size(time_ensemble, 2)/10]);
colormap(hot); 
xlabel('Time (s)','FontSize',14);
colorbar;
title(['Model Ensemble (Last ' num2str(end_count) ' Models)'],'FontSize',14);

subplot(2,3,4);
hist(SigSave(end - end_count + 1:end));
title('Sigma','FontSize',14);

subplot(2,3,5);
hist(NC1Save(end - end_count + 1:end));
title('NoiseCorr1','FontSize',14);

subplot(2,3,6);
hist(NC2Save(end - end_count + 1:end));
title('NoiseCorr2','FontSize',14);

figure(3);
clf;
subplot(2,1,1);
plot(time, P, 'LineWidth', 2);
title('Parent (Reference Phase)');
subplot(2,1,2);
plot(time, D, 'LineWidth', 2);
title('Daughter (Observation)');
xlabel('Time (s)');
sgtitle('Input Data');

%% Display the final model
Dt = time(2) - time(1);
disp('Final Model (Not ensemble):');
fprintf('Location:  %4.2f s; %4.2f s\n', LocSave(end, 1) * Dt, LocSave(end, 2) * Dt);
fprintf('Amplitude: %4.2f;   %4.2f\n', AmpSave(end, 1), AmpSave(end, 2));
fprintf('Width:     %4.2f s; %4.2f s\n', WidSave(end, 1) * Dt, WidSave(end, 2) * Dt);
