% Nikola Janjusevic
% Turbo Coding
close all, clear all;

PLOT = 0;
filename = 'turbo.mat'

M = 2; % constellation order
num_run = 1; % number of runs
num_sym = 15; N = num_sym; % number of symbols
snr_vec = 5; % SNR points

n = 3;
k = 1;
m = 4;
rate = k/n;
% interleave stuff
R=3;C=5;
trellis = poly2trellis(m,[13,15,17],13);
% decoding iterations
l=4;

% SIMULATION
ber_vec = zeros(size(snr_vec));
for ii=1:num_run
    % generate new data each run & pad
    msg = randi([0 1], num_sym, 1);
    % encode
    x1 = turbo_enc(msg,trellis,m); x1 = reshape(x1,n,[]);
    % --- interleave ---
    int_msg = rc_interleave(msg,R,C);
    x2 = turbo_enc(int_msg',trellis,m); x2 = reshape(x2,n,[]);
    x = [x1;x2(2:end,:)];
    % map to constellation (BPSK)
    tx = -2*x+1;
    v = sqrt(1/2)*(randn(size(tx))+1j*randn(size(tx))); % cplx rnd noise
  
    % performs experiment at each SNR
    for jj=1:length(snr_vec)
        rx = tx + 10^(-snr_vec(jj)/20)*v;
        y = reshape(rx,5,[]); ys=y(1,1:N); yp1=y(2,:); yp2=y(4,:);
        int_ys = rc_interleave(ys,R,C);
        Lc = 2*rate*10^(snr_vec(1)/20); Le21 = zeros(N,1);
        for iter=1:l
            % --- D1 ---
            Le12 = mapdec(ys,yp1,trellis,Lc,rc_deinterleave(Le21,R,C),N);
            % --- D2 ---
            Le21 = mapdec(int_ys,yp2,trellis,Lc,rc_interleave(Le12,R,C),N);
        end
        L1 = Lc*ys + rc_deinterleave(Le21,R,C) + Le12';
        msg_hat = L1<0;
        ber_vec(jj) = ber_vec(jj) + biterr(msg,msg_hat')/length(msg);
        fprintf(".");
    end
    fprintf(",\n");
end
ber_vec = ber_vec/num_run
rate = 1/2;
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
