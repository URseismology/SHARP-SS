function D_model = forward_step_transdimensional(P, D, model, prior)
% The forward step of D = P * G in the frequency domain. Returns D in the
% time domain.

% G in time domain
G_model = create_G_from_model_transdimensional(model, prior);

% parameter setup
npts_fft = 2^(1 + nextpow2(length(D)));
fftP = fft(P, npts_fft);

% G in frequency domain
G = [zeros(1,ceil((npts_fft-length(G_model))/2)) ...
    G_model ...
    zeros(1,floor((npts_fft-length(G_model))/2))];

Dtmp = fftshift(ifft(fftP.*fft(G',npts_fft)));
Dtmp = Dtmp';

D_model = Dtmp(1:length(P));