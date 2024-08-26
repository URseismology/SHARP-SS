% reference SS from ricker
fdom = .05; tlen = 60; dt = 0.1;
time = (0:dt:tlen);

[wr,~] = ricker(dt,fdom,2*tlen);  % this generates a ricker wave with a defined dominant frequency
SS_pulse = wr./max(wr);
[SS_ref, ~] = stat(SS_pulse, time, 0);

plot(time, SS_ref);