function [LikeProb, D_model, R_LT ,R_UT ,R_P, LogDetR] = ...
    calculate_like_prob_transdimensional(P, D, model, prior, CdInv_opt, R_LT ,R_UT ,R_P, LogDetR)
% Calculate likelihood probability; also outputs R_LT ,R_UT ,R_P, LogDetR for CdInv
% calculation
%

opts_LT.LT = true;
opts_UT.UT = true;

% create D_model in the frequency domain
D_model = forward_step_transdimensional(P, D, model, prior);

if prior.align

    % align D_model to its peak amplitude (max shift 10 s)
    shifted = 1;
    cut_len = ceil((length(D_model) - 1) * prior.dt / 2 - prior.align); % redefine cut length to 10 s shorter (one side)
    [~, maxind] = max(D_model);
    if abs((maxind - 1) * prior.dt - (length(D_model) - 1) * prior.dt / 2) >= prior.align % don't re-align if > 10 s
        maxind = ceil(length(D_model) / 2);
        shifted = 0;
    end
    midind = ceil(length(D) / 2); % cut from center for D
    D_model = D_model(maxind - round(cut_len/prior.dt):maxind + round(cut_len/prior.dt));
    D = D(midind - round(cut_len/prior.dt):midind + round(cut_len/prior.dt));

    % % debug
    % figure(20);
    % clf;
    % plot(D, 'k-', 'displayname', 'D'); hold on; 
    % plot(D_model, 'r-', 'displayname', 'D_model'); legend;
    % title([num2str(maxind) '  ' num2str(shifted)]);
    % % pause(0.1);

end

if CdInv_opt

    % create CdInv Matrix
    Rrow = zeros(length(D), 1);
    % fit negative part (pre-SS) only
    if prior.negOnly
        Rrow = Rrow(1:round(length(Rrow) / 2));
    end

    for i=1:length(Rrow)

        % The constant in the second cosine is a measure of how many cycles it
        % takes for the decaying-exponential*cosine to decay
        % replace 1.4 by model.NoiseCorr2
        Rrow(i) = exp(-(model.NoiseCorr * (i-1) * prior.dt))...
            * cos(model.NoiseCorr2 * pi * model.NoiseCorr * (i-1) * prior.dt);

    end

    R = toeplitz(Rrow);

    LogDetR = 2 * sum(log(diag(chol(R))));
    [R_LT ,R_UT ,R_P] = lu(R);

end

% calculate difference (P * G - D)
Diff = D_model' - D;
% fit negative part (pre-SS) only
if prior.negOnly
    Diff = Diff(1:round(length(Diff) / 2));
end

Y = linsolve(R_UT, linsolve(R_LT, R_P * Diff, opts_LT), opts_UT);
MahalDist = (Diff.' * Y) / (model.Sigma.^2);
LogCdDeterm = length(D) * log(model.Sigma) + 0.5 * LogDetR;

LikeProb = - LogCdDeterm - MahalDist / 2;

