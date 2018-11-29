% Nikola Janjusevic
% Viterbi Decoding, Hard and Soft
close all, clear all;

PLOT = 0;

M = 2; % constellation order
num_run = 25; % number of runs
num_sym = 1e3; % number of symbols
snr_vec = -5:1:5; % SNR points

n = 2;
k = 1;
m = 5;
rate = k/n;

trellis = poly2trellis(m,[31,27],31);
msg = randi([0 1], 10, 1);
[y final_state] = convenc(msg,trellis);
tail=y(end-1:end);

% forcing the RSCC to the zero state
% mask is the feedback path connections
mask = de2bi(oct2dec(31),m); mask = mask(1:end-1);
for ii=1:(m-1)
    % trellis state in binary
    bstate = de2bi(final_state,m-1);
    % xor'd feedback path in binary, fed into encoder input
    in = mod(sum(mask.*bstate),2);
    [tail,final_state] = convenc(in,trellis,final_state);
    y = [y; tail'];
end
% total of 2*(m-1) tailbits appended
final_state

N = num_symbols;
num_states = trellis.numStates; % NOT related to code rate k/n
num_in = trellis.numInputSymbols; % NOT related to code rate k/n
num_out= trellis.numOutputSymbols;
% map decoding
alpha = nan(num_states,1); alpha(0) = 1;
beta = nan(num_states,1); beta(0) = 0;
gamma = nan(num_states,num_states);
Le = zeros(N,1);

y = reshape(y,2,[]);
for k=1:N
   yk = y(:,k); yks = yk(1); ykp = yk(2:end);
   0
end




%%

ber_vec = zeros(size(snr_vec)); % bit error rate vector
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
