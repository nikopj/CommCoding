% Nikola Janjusevic
% Viterbi Decoding, Hard and Soft
close all, clear all;

PLOT = 0;

M = 2; % constellation order
num_run = 25; % number of runs
num_sym = 1e3; % number of symbols
snr_vec = -5:1:5; % SNR points

for rate = [2/3]
% number of bits for quantization, 1bit is hard decoding
for num_bits = [1 2 3]
ber_vec = zeros(size(snr_vec)); % bit error rate vector
[N,D] = rat(rate);
filename = 'vit_r'+string(N)+string(D)+'_b'+string(num_bits)+'.mat'

% pre-chosen trellis from Wicker
if rate==1/3
    k=1; n=3; K=4;
    trellis = poly2trellis(K,[13 15 17]);
elseif rate==1/2
    k=1; n=2; K=7;
    trellis = poly2trellis(K,[133 171]);
elseif rate==2/3
    k=2; n=3; K=[6 6];
    trellis = poly2trellis(K,[31 46 63; 32 65 61]);
else
    disp('error');
end

% SIMULATION
for ii=1:num_run
  % generate new data each run
  x = randi([0,1],1,num_sym);
  % pad and encode data, force to zero state
  pad = zeros(1, mod(length(x)+max(K)+1,k));
  x = [x pad zeros(1,max(K)+1)];
  % encode
  x_enc = convenc(x, trellis, 0);
  % map to constellation (BPSK)
  tx = -2*x_enc+1;
  v = sqrt(1/2)*(randn(size(tx))+1j*randn(size(tx))); % cplx rnd noise
  
  % performs experiment at each SNR
  for jj=1:length(snr_vec)
    rx = tx + 10^(-snr_vec(jj)/20)*v;
    y = real(rx);
    % cheating by using actual SNR for soft-dec pdfs,
    % could estimate SNR from data, not the purpose of this excercise
    x_hat = myvitdec_soft(y,trellis,num_bits,snr_vec(jj),1,1);
    ber_vec(jj) = ber_vec(jj) + biterr(x,x_hat')/num_sym;
    fprintf(".");
  end
  fprintf(",\n");
end
ber_vec = ber_vec/num_run;

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

save(filename,'-ascii','-double','ber_vec','snr_vec')
end
end