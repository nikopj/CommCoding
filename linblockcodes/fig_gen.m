% figure plotting script
clear all;
load('ham_data')
ber_ham = ber_coded;
ebno_ham = snr_vec + 10*log10(4/7);
load('rm_data')
ber_rm = ber_coded;
ebno_rm = snr_vec + 10*log10(.5);;
load('golay_data')
ber_gol = ber_coded;
ebno_gol = snr_vec + 10*log10(.5);;

figure
semilogy(ebno_ham, ber_ham);
hold on
semilogy(ebno_rm, ber_rm);
hold on
semilogy(ebno_gol, ber_gol);
hold on
ebno = linspace(-2, 12, 100);
ber_theory = qfunc( sqrt( 2*10.^(ebno/10) ) );
semilogy(ebno, ber_theory)

grid on