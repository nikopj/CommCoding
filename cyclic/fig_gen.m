% fig_gen

% ber_5_1: m = 5, t = 1
% n = 31, k=26, rate = .83871
load('bch_5_1');

% ber_5_2: m = 5, t = 2
% n = 31, k=21, rate = .6774
load('bch_5_2');

% ber_5_3: m = 5, t = 3
% n = 31, k=16, rate = .51613
load('bch_5_3');

figure
semilogy(bch_5_1(2,:),bch_5_1(1,:), '-bo', 'MarkerSize', 7, ...
    'LineWidth',2)
hold on
semilogy(bch_5_2(2,:),bch_5_2(1,:), '-rs', 'MarkerSize', 7, ...
    'LineWidth',2)
hold on
semilogy(bch_5_3(2,:),bch_5_3(1,:), '-m+', 'MarkerSize', 7, ...
    'LineWidth',2)
grid on

ebno = linspace(-5,10,100);
ber_theory = qfunc( sqrt( 2*10.^(ebno/10) ) );
hold on
semilogy(ebno,ber_theory, 'k',...
    'LineWidth',2)

legend('t=1, rate=36/31', 't=2, rate=21/31 ', 't=3, rate=16/31', ...
    'uncoded bpsk')

title('BPSK BCH Codes over AWGN, m=5')
xlim([-4, 8])
ylim([.5e-5 .2])
xlabel('EbNo (dB)')
ylabel('BER')