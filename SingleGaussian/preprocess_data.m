function [P, D] = preprocess_data(P, D, lf, hf, prior)


nyq = (0.5 / prior.dt);
fn  = [lf/nyq hf/nyq];
[b,a] = butter(2,fn);

P = double(P);
D = double(D);

b = double(b);
a = double(a);

P = filtfilt(b, a, P);
D = filtfilt(b, a, D);