% fig gen
close all; clear all;

rates = ["13","12","23"];
num_bits = ["1","2","3"];

color = ["r", "b", "m"];
marker = ["+", "o", "*"];
ms=7; lw=1.5;

figure
% uncoded reference
snr = linspace(-5,5,100);
ber_theory = qfunc( sqrt( 2*10.^(snr/10) ) );
semilogy(snr,ber_theory,'k', 'MarkerSize', ms, 'LineWidth', lw);
leg = ['Uncoded BPSK'];
hold on
for ii=1:length(rates)
    r = rates(ii);
    k = r{1}(1); n = r{1}(2);
    rr = str2num(k)/str2num(n);
    
    for jj=1:length(num_bits)
        filename = 'vit_r'+rates(ii)+'_b'+num_bits(jj);
        M = load(filename,'-ascii');
        ber_vec = M(1,:);
        snr_vec = M(2,:);
        
        semilogy(snr_vec+rr,ber_vec,'-'+color(ii)+marker(jj),...
            'MarkerSize',ms,'LineWidth',lw);
        
        leg = [leg 'r='+string(k)+'/'+string(n)+', q='+num_bits(jj)];
        hold on
    end
    
    coded_limit = 10*log10((2^rr -1)/rr);
    plot([coded_limit coded_limit], [1e-6 1e-3], '--'+color(ii))
    leg = [leg 'Shan. Lim, r='+string(k)+'/'+string(n)];
    hold on
end

abs_limit = 10*log10( log(2) );
plot([abs_limit abs_limit], [1e-6 1e-3],'-k')
leg = [leg 'Abs. Shan. Lim'];

xlabel('Eb/N0 (dB)');
xlim([-5 5]);
ylabel('BER');
legend(leg);
grid on
