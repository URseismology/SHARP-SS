data = load('../DataMats/SS_matchfiltered_NoMelt.mat');

P = data.S_seg_corr_avg(1:end-1);
D = data.SS_seg_corr_avg(1:end-1);

time = data.seg_time(1:end-1);

N = length(time);
Dt = time(2) - time(1);
fmax = 1 / (2.0 * Dt);
df = fmax / (N/2);
dw = 2 * pi * df;
w  = dw * [0:N/2,-N/2+1:-1]';

Q0 = 1000;
c0 = 30;
L  = 115000;

qinv0 = 1/Q0;
qinv_w = qinv0;
tstar = L.*qinv_w./c0;

Dwt = attenuation_operator(Q0, c0, L, w, 0, 'mph');
Dwt_td = ifft(Dwt);

fft_P = fft(P);
fft_P_attenuated = fft_P .* Dwt;

P_attenuated = ifft(fft_P_attenuated);
P_attenuated = real(P_attenuated);
P_attenuated = P_attenuated / max(abs(P_attenuated));

figure(1);
clf;
subplot(2,1,1);
plot(time, real(Dwt_td), 'k-');
title(sprintf('Attenuation Operator: t* = %3.2f s', tstar));
subplot(2,1,2);
plot(time, P, 'k-', 'LineWidth', 2, 'DisplayName', 'Original P');
hold on;
% plot(time, D, 'b-', 'LineWidth', 2, 'DisplayName', 'Original D');
% hold on;
plot(time, real(P_attenuated), 'r--', 'LineWidth', 2, 'DisplayName', 'Attenuated P');
legend;

% save to mat file
% save('../DataMats/SS_matchfiltered_NoMelt.mat','P_attenuated','-append');