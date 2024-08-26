% Perform D = P * G to generate synthetics.

% % load P (reference phase) and time
% load(['/Users/evan/Library/CloudStorage/GoogleDrive-evan.z.0920@gmail.com/' ...
%     'My Drive/Research/SS_THBD/SynMats/ricker_P.mat']);
% time = time - 60;

tmp = load(['/Users/evanzhang/Library/CloudStorage/GoogleDrive-evan.z.0920@gmail.com/' ...
    'My Drive/Research/SS_THBD/DataMats/OceanStack.mat']);
time = tmp.time_OS;
P = tmp.SS_mean_10Hz';

% build prior
prior = define_prior(time, P);

% build (customized) model and calculate G_model
TD = 1; % trandimensional (moho+LAB) or not
if ~TD
    loc = 5.0;
    amp = -0.3;
    wid = 1.0;
    model = create_initial_model_customize(prior, loc, amp, wid);
    G_model = create_G_from_model(model, prior);
else
    loc = [5.0 30.0];
    amp = [-0.3 0.04];
    wid = [1.0 3.0];
    model = create_initial_model_TD_customize(prior, loc, amp, wid);
    G_model = create_G_from_model_transdimensional(model, prior);
end

% parameter setup
npts_fft = 2^(1 + nextpow2(length(P)));
fftP = fft(P, npts_fft);

% G in frequency domain
G = [zeros(1,ceil((npts_fft-length(G_model))/2)) ...
    G_model ...
    zeros(1,floor((npts_fft-length(G_model))/2))];

Dtmp = fftshift(ifft(fftP.*fft(G',npts_fft)));

D_model = Dtmp(1:length(P));
G_model = G_model(1:length(P));

% visualize
figure(2);
clf;
subplot(2, 1, 1);
plot(time, P, 'b-', 'LineWidth', 2, 'DisplayName', 'P (Reference)');
hold on;
plot(time, G_model,'k-', 'LineWidth', 2, 'DisplayName', 'G (Model)');
xline(loc, 'k--', 'HandleVisibility', 'off');
xline(-loc, 'k--', 'HandleVisibility', 'off');
legend;
subplot(2, 1, 2);
plot(time, movmean(D_model, 25), 'r-', 'LineWidth', 3, 'DisplayName', 'D = P * G');
xline(loc, 'k--', 'HandleVisibility', 'off');
xline(-loc, 'k--', 'HandleVisibility', 'off');
legend;

% save to file
filename = 'gauss2_OS';
savename = ['/Users/evanzhang/Library/CloudStorage/' ...
    'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/' ...
    'SynMats/' filename '.mat'];
save(savename, 'time', 'P', 'D_model', 'G_model', 'loc', 'amp', 'wid');