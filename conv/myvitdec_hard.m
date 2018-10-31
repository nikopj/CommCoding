function x_hat = myvitdec_hard(y,trellis,init_state,tbb_length)
% y, encoded bits 
% trellis, poly2trellis object
% init_state, number in [0 trellis.numStates-1]
% tbb_length, trace back buffer length

n = trellis.numStates;
k = trellis.numInputSymbols; % NOT related to code rate k/n

states  = trellis.nextStates;
outputs = trellis.outputs;

tb_buffer = zeros(n,tbb_length);
br_metric = zeros(n,1);

y = bi2de(reshape(y,log2(trellis.numSymbols),[])','left-msb');
cur_states = init_state

for i=1:length(y)
    for j=1:length(cur_states)
        br_metric(cur_states(j)) = br_metric_c(cur_states(j)) ...
            + sum(de2bi(bitxor(
    end
    cur_states = 
end

end