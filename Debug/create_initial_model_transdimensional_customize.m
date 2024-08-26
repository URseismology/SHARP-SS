function model = create_initial_model_transdimensional_customize(prior)
% model contains:
% location, amplitude, width: 2-element vectors
% Sigma, NoiseCorr, NoiseCorr2: scalars
%
% prior contains:
% AmpChange, WidthChange, LocChange, SigChange, CorrChange, CorrChange_2,
% tlen (half), dt, std_P,
% minWid, maxWid,

% gaussian parameters
model.location = [round(rand * prior.tlen) round(rand * prior.tlen)];
model.amplitude = [rand rand];
model.width = [round(rand * (prior.maxWid - prior.minWid)) + prior.minWid round(rand * (prior.maxWid - prior.minWid)) + prior.minWid];

% make sure the Gaussian does not fall off the edges
[~, g1_ind] = min(model.location); % find the first gaussian (returns 1 if only one because the other one is inf)
while model.location(1) < (model.width(1) + 1) / 2 || ...
        model.location(1) > prior.tlen - (model.width(1) + 1) / 2 || ...
        model.location(2) < (model.width(2) + 1) / 2 || ...
        model.location(2) > prior.tlen - (model.width(2) + 1) / 2 || ...
        (~isinf(model.location(2)) ...
        && ((model.location(g1_ind) + (model.width(g1_ind) + 1) / 2) ...
        <= (model.location(3 - g1_ind) - (model.width(3 - g1_ind) - 1) / 2)))
    model.location = [round(rand * prior.tlen) round(rand * prior.tlen)];
    model.amplitude = [rand rand];
    model.width = [round(rand * (prior.maxWid - prior.minWid)) + prior.minWid round(rand * (prior.maxWid - prior.minWid)) + prior.minWid];
end

% noise parameters
model.Sigma = prior.std_P;
model.NoiseCorr  = 0.25;
model.NoiseCorr2 = 1.40;