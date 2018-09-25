% Nikola Janjusevic
% Boiler Plate Wireless Channel
% BPSK, AWGN monte-carlo simulation
% plots BER vs. Eb/No with shannon theoretical blah
close all

M = 2; % constellation order
num_run = 10; % number of runs
num_sym = 5e5; % number of symbols
snr_vec = -2:10; % SNR points
ber_vec = zeros(size(snr_vec)); % bit error rate vector

% SIMULATION
for ii=1:num_run
  % generate new data each run
  x = randi([0,1],num_sym,1);
  tx = -2*x+1;
  v = sqrt(1/2)*(randn(num_sym,1)+1j*randn(num_sym,1)); % cplx rnd noise
  
  % performs experiment at each SNR
  for jj=1:length(snr_vec)
    rx = tx + 10^(-snr_vec(jj)/20)*v;
    y  = real(rx)<0;
    ber_vec(jj) = ber_vec(jj) + sum( abs( y-x ) )/num_sym;
  end
end
ber_vec = ber_vec/num_run;

% shannon theoretical ber
rate = 1/2;
ebno = snr_vec - 10*log10(log2(M));
ber_theory_coded = qfunc(sqrt(2*10.^(ebno/10)/rate));

% shannon limit
abs_limit = 10*log10( log(2) );
coded_limit = 10*log10( (2^rate -1)/rate);

% plotting
figure
% bers
semilogy(ebno, ber_vec, '-bo', ebno, ber_theory_coded, 'r')
% limits
ymin = min(ber_vec(ber_vec~=0));
line([abs_limit abs_limit], [ymin 1e-3], ...
    'color', 'black', 'linestyle', '--')
line([coded_limit coded_limit], [ymin 1e-3], ...
    'color', [0.5 0 1], 'linestyle', '--')
title('BPSK over AWGN Channel')
xlim([ebno(1) ebno(end)])
ylim([ymin .5])
xlabel('Eb/N0 (dB)')
ylabel('BER')
%ylim([1e-5 1e-1])
legend('simulation, rate: 1','theory, rate: '+string(rate), ...
    'abs. Shannon limit', 'Shannon limit, rate: '+string(rate))
grid on