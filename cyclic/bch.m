% Nikola Janjusevic
% BCH, BPSK, AWGN monte-carlo simulation
% plots BER vs. Eb/No with shannon theoretical blah
close all; clear all;

M = 2; % constellation order
num_run = 10; % number of runs
num_sym = 1e4; % number of symbols
snr_vec = -2:8; % SNR points
ber_coded = zeros(size(snr_vec)); % bit error rate vector

% BCH defined by m, t
m = 3;
t = 1;
n = 2^m-1;

[G, H] = myBCHgen(m,t);
k = size(G,1);
rate = k/n;

% to save ber and ebno
filename = 'bch_'+string(m)+'_'+string(t)

pad = zeros(mod(k-num_sym, k), 1);
% SIMULATION
for ii=1:num_run
  % generate new data each run
  x = [randi([0,1],num_sym,1); pad];
  % reshapes message data into chunks of len k, then codes them with G
  % and reshapes them back into one column.
  x_enc = mod( reshape( G'*reshape( x, k, [] ) , [], 1), 2 );
  % tx is mapped transmission ( ie -1,1)
  tx = -2*x_enc+1;
  % cplx wg noise
  v = sqrt(1/2)*( randn(length(tx), 1) + 1j*randn(length(tx), 1) ); 
  
  % performs experiment at each SNR
  for jj=1:length(snr_vec)
    rx = tx + 15^(-snr_vec(jj)/20)*v;
    % each coded word in columns
    y_enc  = reshape(real(rx)<0, n, []);
    % need to reverse for right order in galois field polynomial
    ygf_enc = gf( flipud(y_enc), m ); 
    % syndrom of each message in columns
    S = mod(H*y_enc, 2);
    E = myBCHerror(S, ygf_enc, m, t);
    y_hat = mod(y_enc + E, 2);
    x_hat = reshape( [eye(k) zeros(k,n-k)]*y_hat, [], 1 );
    
    ber_coded(jj) = ber_coded(jj) + sum( abs( x_hat-x ) )/length(x);
    fprintf(".");
  end
  fprintf(",\n");
end
ber_coded = ber_coded/num_run

% theoretical ber
ebno = snr_vec - 10*log10(log2(M)); %+ 10*log10(rate);
ber_theory = qfunc( sqrt( 2*10.^(ebno/10) ) );

% shannon limit
abs_limit = 10*log10( log(2) );
coded_limit = 10*log10( (2^rate -1)/rate);

% plotting
figure
% bers
semilogy(ebno, ber_coded, '-bo', ... 
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