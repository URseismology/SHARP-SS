function prior = define_prior(time, P)

% prior contains:
% AmpChange, WidthChange, LocChange, SigChange, CorrChange, CorrChange_2,
% tlen (half), dt, std_P,
% minWid, maxWid,


prior.WidthChange  = 0.2;
prior.LocChange    = 0.2;
prior.AmpChange    = 0.002;
prior.SigChange    = 0.0025;
prior.CorrChange   = 0.000125;
prior.CorrChange_2 = 0.012;

prior.tlen  = round((length(time) + 1) / 2); % half length in sample points
prior.dt    = time(2) - time(1);
prior.std_P = std(P);

prior.maxWid = prior.tlen / 4;
prior.minWid = 3; % 3 % in sample points; default 3 (min req. for Gaussian)

prior.minAmp = -0.5;
prior.maxAmp = 0.5;

prior.minLoc  = 1.0; % in seconds

prior.negOnly = 0; % fit precursor part only

prior.FOM     = 0; % force oceanic Moho at 4.0 s with the amplitude of 0.2

prior.align   = 0; % align at peak for each new D; number here is max shift in seconds