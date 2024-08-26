% test a 1-D search for LAB depth in NoMelt data ...

%%%% NoMelt Data (Match Filtered)

data = load(['/Users/evanzhang/Library/CloudStorage/' ...
    'GoogleDrive-evan.z.0920@gmail.com/My Drive/Research/SS_THBD/DataMats/' ...
    'SS_matchfiltered_NoMelt.mat']);
% P = data.S_seg_corr_avg(1:end-1);
P = data.P_attenuated;
D = data.SS_seg_corr_avg(1:end-1);
time = data.seg_time(1:end-1);

%%%% Define model and prior

prior = define_prior(time, P);

model.location   = [40 inf];
model.amplitude  = [-0.3 0.02];
model.width      = [10 30];
model.Sigma      = std(P);
model.NoiseCorr  = 0.25;
model.NoiseCorr2 = 1.40;

%%%% Loop through different LAB location

LAB_loc = (100:1:450);
LP_vec  = zeros(1, length(LAB_loc));

for i = 1:length(LAB_loc)

    % update model
    model.location(2) = LAB_loc(i);

    % calculate Cd_Inv at step 1
    if i == 1
        [model.LikeProb, D_model, R_LT ,R_UT ,R_P, LogDetR] = ...
            calculate_like_prob_transdimensional(P, D, model, prior, 1);
    else
        [model.LikeProb, D_model, ~, ~, ~, ~] = ...
            calculate_like_prob_transdimensional(P, D, model, prior, 0, R_LT ,R_UT ,R_P, LogDetR);
    end

    LP_vec(i) = model.LikeProb;

end

figure(1);
clf;
% plot(LAB_loc * prior.dt *4.1 / 2, LP_vec, 'k-', 'LineWidth', 2);
% xlabel("LAB Depth (km)");
plot(LAB_loc * prior.dt, LP_vec, 'k-', 'LineWidth', 2);
xlabel("LAB Location (s)");
ylabel("Likelihood Probability");
title(['Align: ' num2str(prior.align) '; NegOnly: ' num2str(prior.negOnly)]);