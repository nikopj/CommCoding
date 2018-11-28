function x_hat = mymap_dec(y,trellis,num_bits,SNR,init_state,final_state)
% y, encoded bits 
% trellis, poly2trellis object
% init_state, number in [0 trellis.numStates-1]
% tbb_length, trace back buffer length

% ASSUMING BPSK @ +/-1!!!!

num_st = trellis.numStates; % NOT related to code rate k/n
num_in = trellis.numInputSymbols; % NOT related to code rate k/n
num_out= trellis.numOutputSymbols;

states  = trellis.nextStates;
bin_outs = nan(size(trellis.outputs));
bin_outs = repmat(bin_outs,1,1,log2(num_out));
for ii=1:size(bin_outs,1)
    for jj=1:size(bin_outs,2)
        bin_outs(ii,jj,:) = de2bi(trellis.outputs(ii,jj),log2(num_out),...
            'left-msb');
    end
end

% starting at init_state, so all other branch metrics 
% are weighed infinitely negatively
br_metric = -inf(num_st,1);
br_metric(init_state) = 0;

y = reshape(y,log2(num_out),[]);
tb_buffer = nan(num_st,size(y,2));

% soft_decoding stuff
range = 0.8*[min(min(y)) max(max(y))];
M = nan(2,2^num_bits);
r = linspace(range(1),range(2),2^num_bits+1);
sig = 10^(-SNR/20);
for ii=[1 2]
    % mean of bpsk pdf
    mu = -2*(ii-1)+1;
    for jj=1:(2^num_bits)
        % integrating probability in small section
        % Q(min) - Q(max)
        M(ii,jj) = log(qfunc((r(jj)-mu)/sig) - qfunc((r(jj+1)-mu)/sig));
    end
end

% quantize bits
y_quant = nan(size(y));
for ii=1:size(y,2)
    y_quant(:,ii) = sum(y(:,ii)>r(1:end-1),2);
end
y_quant(y_quant==0) = 1;

for ii=1:size(y,2)
    br_temp = -inf(size(br_metric));
    tb_temp = nan(size(tb_buffer));
    for jj=1:num_st
        % find all possible conditions for entering the current state,
        % jj, where I are the previous states, and J are the inputs
        [I,J] = find( states==(jj-1) );
        metrics = -inf(num_in,1);
        % for number of entering branches
        for kk=1:length(I)
            out = reshape(bin_outs(I(kk),J(kk),:), [],1);
            % lookup-table of conditional metrics in M, 
            % take diagonal for matlab reasons,
            % sum the metrics for the word
            metrics(kk) = sum(diag(M(out+1,y_quant(:,ii))));
        end       
        % find survivor branch
        [~,K] = max(br_metric(I) + metrics);
        % add to branch metric
        br_temp(jj) = br_metric(I(K)) + metrics(K);
        % add input to traceback buffer
        tb_temp(jj,:) = tb_buffer(I(K),:); tb_temp(jj,ii) = J(K)-1;
    end
    br_metric = br_temp;
    tb_buffer = tb_temp;
end

x_hat = reshape(de2bi(tb_buffer(final_state,:), ... 
            log2(num_in), 'left-msb')', [], 1);
end