function model = create_initial_model_TD_customize(prior, loc, amp, wid)
% model contains:
% location, amplitude, width: 2-element vectors
% Sigma, NoiseCorr, NoiseCorr2: scalars
%
% prior contains:
% AmpChange, WidthChange, LocChange, SigChange, CorrChange, CorrChange_2,
% tlen (half), dt, std_P,
% minWid, maxWid,

% gaussian parameters
model.location = round(loc / prior.dt);
model.amplitude = amp;
model.width = round(wid / prior.dt);

% noise parameters
model.Sigma = prior.std_P;
model.NoiseCorr  = 0.25;
model.NoiseCorr2 = 1.40;