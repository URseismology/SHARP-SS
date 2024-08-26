% Add realistic noise from noise paramaterization in Kolb and Lekic (2014)

% load P (reference phase) and time
tmp = load(['/Users/evan/Library/CloudStorage/GoogleDrive-evan.z.0920@gmail.com/' ...
    'My Drive/Research/SS_THBD/SynMats/ricker_P.mat']);
time = tmp.time - 60;
P = tmp.Z;

% build prior
prior = define_prior(time, P);
data_ct = 1000;

%% Generate noise

% Define noise paramters and construct covariance matrix
lambda = 0.2;
omega0 = 4.4;

lenT = length(time);

CovMat = zeros(lenT);

for i = 1:lenT
    for j = 1:lenT
        CovMat(i, j) = exp(-lambda * abs(time(j) - time(i))) ...
            * cos(lambda * omega0 * abs(time(j) - time(i)));
    end
end

% Generate noise
L = chol(CovMat, 'lower');
noise = L * randn(lenT, data_ct);
noise = (normalize(noise, 1, "range") - 0.5) * 2; % normalize

%% Genrate G (and D)

allX = zeros(data_ct, length(P), 2);
allY = zeros(data_ct, length(P));
allY_params = zeros(data_ct, 6); % 3 - loc, amp, wid * 2

for i = 1:data_ct

    % build (customized) model and calculate G_model

    TD = round(rand); % trandimensional (moho+LAB) or not - randomly

    if ~TD
        loc = rand*8+12.0; % 12-20 randomly
        amp = rand*0.25+0.25; % 0.25-0.5 randomly
        wid = rand+0.5; % 0.5-1.5 randomly
        model = create_initial_model_customize(prior, loc, amp, wid);
        G_model = create_G_from_model(model, prior);
    else
        loc = [rand*5+2.0 rand*8+12.0]; % 2-7 and 12-20 randomly
        amp = [-(rand*0.25+0.25) (rand*0.15+0.05)]; % -(0.25-0.5) and 0.05-0.2 randomly
        wid = [rand+0.5 rand*2+2]; % 0.5-1.5 and 2-4 randomly
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

    thisNoise = (rand * 0.5) * noise(:, i); % noise level chosen randomly 0-0.5
    D_model = D_model + thisNoise;

    allX(i, :, 1) = P;
    allX(i, :, 2) = D_model;
    allY(i, :) = G_model;

    if TD
        allY_params(i, :) = [loc(1); amp(1); wid(1); loc(2); amp(2); wid(2)];
    else
        allY_params(i, :) = [loc(1); amp(1); wid(1); 0; 0; 0];
    end

    % visualize
    figure(2);
    clf;
    subplot(2, 1, 1);
    plot(time, P, 'b-', 'LineWidth', 2, 'DisplayName', 'P (Reference)');
    hold on;
    plot(time, G_model,'k-', 'LineWidth', 2, 'DisplayName', 'G (Model)');
    legend;
    subplot(2, 1, 2);
    plot(time, D_model, 'r-', 'LineWidth', 3, 'DisplayName', 'D = P * G');
    legend;

end



% Example code to save the input data X to a CSV file using writematrix
% Reshape X to 2D matrix for CSV
allX_reshaped = reshape(allX, [], size(allX, 3));

% Save data to CSV file
savedir = ['/Users/evan/Library/CloudStorage/GoogleDrive-evan.z.0920@gmail.com/' ...
    'My Drive/Research/SS_THBD/NN_test/SynData/Test_02'];

writematrix(allX_reshaped, [savedir 'allX_reshaped.csv']);
writematrix(allY, [savedir 'allY.csv']);
writematrix(allY_params, [savedir 'allY_params.csv']);
