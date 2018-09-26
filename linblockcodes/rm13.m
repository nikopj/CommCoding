% Nikola Janjusevic
% Boiler Plate Wireless Channel
% BPSK, AWGN monte-carlo simulation
% plots BER vs. Eb/No with shannon theoretical blah
close all; clear all;

M = 2; % constellation order, BPSK
num_run = 10; % number of runs
num_sym = 5e5; % number of symbols
snr_vec = 1:20; % SNR points
ber_coded   = zeros(size(snr_vec));

global m
m = 5;
r = 1;
n = 2^m;
global k
k = 1+m;
H = hadamard(2);
global Hm
Hm = eye(n);
for j=1:m
   Hm = Hm * kron( kron( eye(2^(m-j)), H), eye(2^(j-1)) );
end
G = Grm(r,m);
rate = k/n;

pad = zeros(k - mod(num_sym, k), 1);

% SIMULATION
for ii=1:num_run
  % generate new data each run
  x = [randi([0,1],num_sym,1); pad];
  % reshapes message data into chunks of 4, then codes them with G
  % and reshapes them back into one column.
  x_enc = mod( reshape( G'*reshape( x, k, [] ) , [], 1), 2 );
  % tx is mapped transmission ( ie -1,1)
  tx = -2*x_enc+1;
  % cplx wg noise
  v = sqrt(1/2)*( randn(length(tx), 1) + 1j*randn(length(tx), 1) ); 
  
  % performs experiment at each SNR
  for jj=1:length(snr_vec)
    rx = tx + (1/rate)*10^(-snr_vec(jj)/20)*v;
    % each coded word in columns
    y_enc  = reshape(real(rx)<0, n, []);
    xx = decode(y_enc);
    x_hat = reshape(xx, [], 1);
    
    ber_coded(jj) = ber_coded(jj) + sum( abs( x_hat-x ) )/length(x);
    fprintf(".");
  end
  fprintf(",\n");
end
ber_coded = ber_coded/num_run

% theoretical ber
ebno = snr_vec - 10*log10(log2(M)) + 10*log10(rate);
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
ylim([ymin 1])
xlabel('Eb/N0 (dB)')
ylabel('BER')
%ylim([1e-5 1e-1])
legend('simulation','theory, uncoded', ...
    'abs. Shannon limit', 'Shannon limit, rate: '+string(rate))
grid on

function decd = decode(encd)
global Hm; global m; global k;
vi0 = zeros(k-1, size(encd,2));
w0 = 2*encd - 1;
wm = Hm'*w0;
i0 = zeros(1, size(encd,2));
find_arg = (abs(wm) == max( abs(wm), [], 1));
for j=1:size(encd,2)
    [i0(1,j), ~] = find( find_arg(:,j), 1, 'first' );
end
s = zeros(1, size(encd,2));
for j=1:size(encd,2)
    s(1,j) = wm(i0(1,j),j) > 0;
end
v = de2bi(i0-1)';
vi0(1:size(v,1),:) = v;
decd = [s ; vi0];
end