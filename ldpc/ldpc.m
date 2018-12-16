% Nikola Janjusevic
% LDPC Decoding. *ldpcCode function provided by Karol
close all, clear all;

PLOT = 1;
filename = 'LDPC.mat';

M = 2; % constellation order
num_run = 100; % number of runs
num_sym = 1e3; % number of symbols
snr_vec = -1:1:12; % SNR points

[H,A,B,T,C,D,E] = ldpcCode(); % *Karol's function
n = size(H,2); k = n-size(H,1);
rate = k/n;

[I,J] = find(H~=0);
MCi = zeros(length(I)); MVj = zeros(length(I));
LCi = nan(1,n); LVj = nan(1,n-k);
for ii=1:n
    c = setC(ii,I,J); LCi(ii) = length(c);
    MCi(ii,1:LCi(ii)) = c;
end
for jj=1:(n-k)
    v =  setV(jj,I,J); LVj(jj) = length(v);
    MVj(jj,1:LVj(jj)) = v;
end

% SIMULATION
ber_vec = zeros(size(snr_vec));
for ii=1:num_run
  % generate new data each run & pad
  x = randi([0,1],num_sym,1); x = [x; zeros(k-mod(num_sym,k),1)];
  % encode data
  x_enc = ldpc_enc(x,H,A,B,T,C,D,E);
  % map to constellation (BPSK)
  tx = -2*x_enc+1;
  v = sqrt(1/2)*(randn(size(tx))+1j*randn(size(tx))); % cplx rnd noise
  
  % performs experiment at each SNR
  for jj=1:length(snr_vec)
    rx = tx + 10^(-snr_vec(jj)/20)*v;
    [M,y] = soft_output(real(rx),1,snr_vec(jj));
    x_hat = ldpc_msgpass2(y,H,M,MCi,LCi,MVj,LVj);
    ber_vec(jj) = ber_vec(jj) + biterr(x,x_hat)/length(x);
    fprintf(".");
  end
  fprintf(",\n");
end
ber_vec = ber_vec/num_run

if(PLOT)
% theoretical ber
ber_theory = qfunc( sqrt( 2*10.^(snr_vec/10) ) );

% shannon limit
abs_limit = 10*log10( log(2) );
coded_limit = 10*log10( (2^rate -1)/rate);

% plotting
figure
% bers
semilogy(snr_vec, ber_vec, '-bo', ... 
    snr_vec, ber_theory, 'r')
% limits
ymin = min(ber_vec(ber_vec~=0));
line([abs_limit abs_limit], [-inf 1e-3], ...
    'color', 'black', 'linestyle', '--')
line([coded_limit coded_limit], [-inf 1e-3], ...
    'color', [0.5 0 1], 'linestyle', '--')
title('BPSK over AWGN Channel')
xlim([snr_vec(1) snr_vec(end)])
%ylim([ymin .5])
xlabel('SNR (dB)')
ylabel('BER')
%ylim([1e-5 1e-1])
legend('simulation','theory, uncoded', ...
    'abs. Shannon limit', 'Shannon limit, rate: '+string(rate))
grid on
end

save(filename,'-mat','-double','ber_vec','snr_vec')