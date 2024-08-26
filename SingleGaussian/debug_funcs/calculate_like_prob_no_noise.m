function [LikeProb, R_LT ,R_UT ,R_P, LogDetR] = ...
    calculate_like_prob_no_noise(P, D, model, prior, CdInv_opt, R_LT ,R_UT ,R_P, LogDetR)
% Calculate likelihood probability; also outputs R_LT ,R_UT ,R_P, LogDetR for CdInv 
% calculation
%

opts_LT.LT = true; 
opts_UT.UT = true; 

if CdInv_opt

    % create CdInv Matrix
    Rrow = zeros(length(D), 1);

    for i=1:length(D)

        % The constant in the second cosine is a measure of how many cycles it
        % takes for the decaying-exponential*cosine to decay
        Rrow(i) = exp(-(model.NoiseCorr * (i-1) * prior.dt))...
            * cos(1.4 * pi * model.NoiseCorr * (i-1) * prior.dt);

    end

    R = toeplitz(Rrow);

    LogDetR = 2 * sum(log(diag(chol(R))));
    [R_LT ,R_UT ,R_P] = lu(R);

end

% create D_model in the frequency domain
D_model = forward_step(P, D, model, prior);

% calculate difference (P * G - D)
Diff = D_model' - D;

MahalDist = (Diff.' * Diff) / (model.Sigma.^2);
LogCdDeterm = length(D) * log(model.Sigma) + 0.5 * 2 * sum(log(ones(length(D), 1)));

LikeProb = -LogCdDeterm - MahalDist / 2;