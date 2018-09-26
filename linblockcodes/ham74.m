% Nikola Janjusevic
% Boiler Plate Wireless Channel
% BPSK, AWGN monte-carlo simulation
% plots BER vs. Eb/No with shannon theoretical blah
close all; clear all

M = 2; % constellation order, BPSK
num_run = 1; % number of runs
num_sym = 5e5; % number of symbols
snr_vec = 1:16; % SNR points
ber_coded   = zeros(size(snr_vec));

m = 3;
global n
n = 2^m-1;
k = n - m;
rate = k/n;

Q = [1 0 1 1; 1 1 1 0; 0 1 1 1];
H = [eye(3) Q];
G = [Q' eye(4)];

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
    % syndrom of each message in columns
    S = mod(H*y_enc, 2);
    E = e(S);
    y_hat = mod(y_enc + E, 2);
    x_hat = reshape( y_hat(k:end,:), [], 1 );
    
    ber_coded(jj) = ber_coded(jj) + sum( abs( x_hat-x ) )/length(x);
    fprintf(".");
  end
  fprintf(",\n");
end
ber_coded = ber_coded/num_run;

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
ylim([ymin .5])
xlabel('Eb/N0 (dB)')
ylabel('BER')
%ylim([1e-5 1e-1])
legend('simulation','theory, uncoded', ...
    'abs. Shannon limit', 'Shannon limit, rate: '+string(rate))
grid on

% error vector generator, able to take in a matrix of syndromes
% where each syndrome is a column of S
function E = e(S)
global n
E = zeros(n, size(S,2));
E(1,:) = (S(1,:)).*(~S(2,:)).*(~S(3,:));
E(2,:) = (~S(1,:)).*(S(2,:)).*(~S(3,:));
E(3,:) = (~S(1,:)).*(~S(2,:)).*(S(3,:));
E(4,:) = (S(1,:)).*(S(2,:)).*(~S(3,:));
E(5,:) = (~S(1,:)).*(S(2,:)).*(S(3,:));
E(6,:) = (S(1,:)).*(S(2,:)).*(S(3,:));
E(7,:) = (S(1,:)).*(~S(2,:)).*(S(3,:));
end