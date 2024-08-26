% debug: test all actions

prior.dt = 0.1;

prior.WidthChange  = 0.2;
prior.LocChange    = 0.5;
prior.AmpChange    = 0.05;
prior.SigChange    = 0.0025;
prior.CorrChange   = 0.000125;
prior.CorrChange_2 = 0.012;

totalSteps = 10000;

% Location change

loc = zeros(1, totalSteps);
loc(1) = 100;

for i = 2:totalSteps
    loc(i) = loc(i-1) + round(prior.LocChange * randn/prior.dt);
end

figure(1);
subplot(3,1,1);
plot(loc);

% Amplitude change

amp = zeros(1, totalSteps);
amp(1) = 0.1;

for i = 2:totalSteps
    amp(i) = amp(i-1) + prior.AmpChange * randn;
    amp(i) = min(max(amp(i), prior.minAmp), prior.maxAmp);
end

figure(1);
subplot(3,1,2);
plot(amp);

% Width change

wid = zeros(1, totalSteps);
wid(1) = 40;

for i = 2:totalSteps
    wid(i) = wid(i-1) + round(prior.WidthChange * randn/prior.dt);
    wid(i) = min(max(wid(i), prior.minWid), prior.maxWid);
end

figure(1);
subplot(3,1,3);
plot(wid);
