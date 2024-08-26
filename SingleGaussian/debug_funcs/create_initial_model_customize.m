function model = create_initial_model_customize(prior, loc, amp, wid)
% model contains:
% location, amplitude, width, Sigma, NoiseCorr, NoiseCorr2
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