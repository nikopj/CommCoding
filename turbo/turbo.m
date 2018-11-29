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
[x final_state] = convenc(msg,trellis);
tail=x(end-1:end);

% forcing the RSCC to the zero state
% mask is the feedback path connections
mask = de2bi(oct2dec(31),m); mask = mask(1:end-1);
for ii=1:(m-1)
    % trellis state in binary
    bstate = de2bi(final_state,m-1);
    % xor'd feedback path in binary, fed into encoder input
    in = mod(sum(mask.*bstate),2);
    [tail,final_state] = convenc(in,trellis,final_state);
    x = [x; tail'];
end
% total of 2*(m-1) tailbits appended
final_state

tx = -2*x+1;
rx = tx;

N = num_sym;
num_states = trellis.numStates;
num_in = trellis.numInputSymbols; 
num_out= trellis.numOutputSymbols;
% map decoding
S = trellis.nextStates; S_minus=S(:,2); S_plus=S(:,1);
% alpha and beta are really alpha,beta tilda
alpha = zeros(num_states,1); alpha(1) = 1;
beta = zeros(num_states,1); beta(1) = 0;
gamma = nan(num_states,num_states);
gamma_e = nan(num_states,num_states);
Le = zeros(N,1);
snr = 20; %db
Lc = 10^(snr/10)*8;

% parity bit generation
x_parity = nan(num_states,num_states);
uk = nan(num_states,num_states);
for ii=1:num_states
    for input=1:num_in
        enc = de2bi(trellis.outputs(ii,input),log2(num_out),'left-msb');
        par = bi2de(enc(2:end));
        jj  = S(ii,input)+1;
        x_parity(ii,jj) = par;
        % input value to the encoder
        uk(ii,jj) = -2*(input-1)+1;
    end
end
x_parity = -2*x_parity+1;
% valid transitions from one state to the next
[s_prime,s]=find(isnan(x_parity)==1);
%%
y = reshape(rx,2,[]);
for k=1:N
   yk = y(:,k); yks = yk(1); ykp = yk(2:end);
   % compute gammas
   for t=1:length(s)
       xkp = x_parity(s_prime(t),s);
       u   = uk(s_prime(t),s(t)); 
       gamma_e(s_prime(t),s(t))   = exp(0.5*Lc*ykp*xkp);
       gamma(s_prime(t),s(t))   = exp(0.5*u*(Le+Lc*yks))*gamma_e(s_prime(t),s(t));
   end
   % compute alphas
   % numerator is single sum over primes, denominator is double
   % sum over primed and unprimed
   den = 0;
   for t=1:length(s)
       num=0;
       for tp=1:length(s_prime)
           num = num + alpha(tp)*gamma(tp,t);
           den = den + alpha(tp)*gamma(tp,t);
       end
       alpha(t) = num;
   end
   alpha = alpha/den;
end
% betas betas betas
for k=N:-1:2
    den = 0;
    for tp=1:length(s_prime)
       num=0;
       for t=1:length(s)
           num = num + beta(t)*gamma(tp,t);
           den = den + alpha(tp)*gamma(tp,t);
       end
       beta(t) = num;
    end
    alpha = alpha/den;
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
    x = real(rx);
    % cheating by using actual SNR for soft-dec pdfs,
    % could estimate SNR from data, not the purpose of this excercise
    x_hat = myvitdec_soft(x,trellis,num_bits,snr_vec(jj),1,1);
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
