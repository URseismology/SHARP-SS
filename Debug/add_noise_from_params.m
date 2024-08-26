% Add realistic noise from noise paramaterization in Kolb and Lekic (2014)

%% Generate noise
% Define directories and noise counts
original_dir = ['/Users/evanzhang/Library/CloudStorage/' ...
    'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/SynMats/'];

output_dir   = ['/Users/evanzhang/Library/CloudStorage/' ...
    'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/SynMats/' ...
    'Noise/'];

modname = 'gauss2_OS';
orig_mat = load([original_dir modname '.mat']);

noise_ct = 50; % 50 sets of random noise (noise parameters are the same)

% Define noise paramters and construct covariance matrix
lambda = 0.2;
omega0 = 4.4;

time = orig_mat.time;
dt   = time(2) - time(1);
lenT = length(time);

G_model = orig_mat.G_model;
P = orig_mat.P;
amp = orig_mat.amp;
loc = orig_mat.loc;
wid = orig_mat.wid;

CovMat = zeros(lenT);

for i = 1:lenT
    for j = 1:lenT
        CovMat(i, j) = exp(-lambda * abs(time(j) - time(i))) ...
            * cos(lambda * omega0 * abs(time(j) - time(i)));
    end
end

% Generate noise
L = chol(CovMat, 'lower');
noise = L * randn(lenT, noise_ct);
noise = (normalize(noise, 1, "range") - 0.5) * 2; % normalize

%% Add noise
% Loop through all synthetics in orignal dir and add noise

if ~exist([output_dir modname], 'dir')
    mkdir([output_dir modname])
end

for inoise = 1:noise_ct

    thisNoise = 0.1 * noise(:, inoise);
    D_model = orig_mat.D_model + thisNoise;

    save([output_dir modname '/' ...
        modname...
        '_Noise' num2str(inoise, '%02d') '.mat'],...
        'D_model', 'G_model', 'P', 'amp', 'loc', 'time', 'wid');

end