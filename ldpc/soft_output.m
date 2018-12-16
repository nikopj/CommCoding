function [M,y_quant] = soft_output(y,num_bits,SNR)
% soft_decoding stuff
range = 0.8*[min(min(y)) max(max(y))];
M = nan(2,2^num_bits);
if num_bits==1
   r =  sum(range)/2;
elseif num_bits==2
    r = [range(1) 0 range(2)];
else
    r1 = linspace(range(1),0,2^(num_bits-1));
    r2 = linspace(0,range(2),2^(num_bits-1));
    r = [r1(1:end-1) 0 r2(2:end)];
end
r = [-inf r inf];
sig = 10^(-SNR/20);
for ii=[1 2]
    % mean of bpsk pdf
    mu = -2*(ii-1)+1;
    for jj=1:(2^num_bits)
        % integrating probability in small section
        % Q(min) - Q(max)
        M(ii,jj) = qfunc((r(jj)-mu)/sig) - qfunc((r(jj+1)-mu)/sig);
    end
    %M(ii,end) = qfunc((r(jj)-mu)/sig);
end

% quantize bits
y_quant = nan(size(y));
for ii=1:size(y,2)
    y_quant(:,ii) = sum(y(:,ii)>r(1:end-1),2);
end
y_quant(y_quant==0) = 1;
y_quant(y_quant==2^num_bits+1) = 2^num_bits; 
end