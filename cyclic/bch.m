% BCH

% Nikola Janjusevic
% Boiler Plate Wireless Channel
% BPSK, AWGN monte-carlo simulation
% plots BER vs. Eb/No with shannon theoretical blah
close all; clear all;

M = 2; % constellation order
num_run = 5; % number of runs
num_sym = 5e5; % number of symbols
snr_vec = 0:20; % SNR points
ber_coded = zeros(size(snr_vec)); % bit error rate vector

m = 4;
t = 2;
n = 2^m-1;

a = gf(2,m);
jj = 1:2:(2*t-1);
phi_vec = sym(zeros(1,length(jj)));
i=1;
for j=jj
    mp = minpol(a^j);
    phi_vec(i) = poly2sym( cast(mp.x,'like',1) );
    i=i+1;
end
lcm_sym = lcm(phi_vec);
lcm_poly = mod(sym2poly(lcm_sym), 2);
g = gf(lcm_poly,m);
k = n - length(g) + 1;

% making generator matrix, G, from circular shifts of generator polynomial
% vector, g
pad = zeros(1, n - mod(length(g), n) );
gx = [pad g.x];
G = zeros(k,n)*nan;
for i=1:k
    gx = circshift(gx,1);
    G(i,:) = gx;
end

H = zeros(length(jj)*m, n)*nan;
i=1;
for j=jj
    b = 0:(n-1);
    H(i:(i+m-1),:) = gftuple((b.*j)',m)';
    i=i+m;
end
msg = gf(randi([0 1],1,k), m);
Gg = gf(G,m);
Hg = gf(H,m);
v = msg*Gg
s = v*Hg'


%%

pad = zeros(k - mod(num_sym, k), 1);

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
    rx = tx + (1/rate)*10^(-snr_vec(jj)/20)*v;
    % each coded word in columns
    y_enc  = reshape(real(rx)<0, n, []);
    % syndrom of each message in columns
    S = mod(H*y_enc, 2);
    E = e(S);
    y_hat = mod(y_enc + E, 2);
    x_hat = reshape( [zeros(k,k) eye(k)]*y_hat, [], 1 );
    
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
disp('shit');
end