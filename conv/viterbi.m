% Nikola Janjusevic
% Boiler Plate Wireless Channel
% BPSK, AWGN monte-carlo simulation
% plots BER vs. Eb/No with shannon theoretical blah
close all, clear all;

M = 2; % constellation order
num_run = 1; % number of runs
num_sym = 1e4; % number of symbols
snr_vec = -5:2:10; % SNR points
ber_vec = zeros(size(snr_vec)); % bit error rate vector

konst = [4 3];
g = [4 5 17;7 4 2];
k = length(konst);
n = size(g,2);
rate = k/n;
trellis = poly2trellis(konst,g);

% SIMULATION
for ii=1:num_run
  % generate new data each run
  x = randi([0,1],1,num_sym);
  
  % pad and encode data
  pad = zeros(1, mod(length(x)+max(konst)-1,k));
  x = [x pad zeros(1,max(konst)-1)];
  x_enc = convenc(x, trellis);
  % map to constellation
  tx = -2*x+1;
  v = sqrt(1/2)*(randn(num_sym,1)+1j*randn(num_sym,1)); % cplx rnd noise
  
  % performs experiment at each SNR
  for jj=1:length(snr_vec)
    rx = tx + 10^(-snr_vec(jj)/20)*v;
    y  = real(rx)<0;
    x_hat = myvitdec_hard(x_enc,trellis,1,1);
    ber_vec(jj) = ber_vec(jj) + biterr(x,x_hat')/num_sym;
    fprintf(".");
  end
  fprintf(",\n");
end
ber_vec = ber_vec/num_run;

% theoretical ber
ebno = snr_vec - 10*log10(log2(M)) - 10*log10(rate);
ber_theory = qfunc( sqrt( 2*10.^(ebno/10) ) );

% shannon limit
abs_limit = 10*log10( log(2) );
coded_limit = 10*log10( (2^rate -1)/rate);

% plotting
figure
% bers
semilogy(ebno, ber_vec, '-bo', ... 
    ebno, ber_theory, 'r')
% limits
ymin = min(ber_coded(ber_coded~=0));
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
legend('simulation','theory, uncoded', ...
    'abs. Shannon limit', 'Shannon limit, rate: '+string(rate))
grid on

save(filename,'-ascii','-double','ber_coded','ebno')