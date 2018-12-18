function x = turbo_enc(msg,trellis,m)
[x,final_state] = convenc(msg,trellis);
tail=x(end-1:end);

% forcing the RSCC to the zero state
% mask is the feedback path connections
mask = de2bi(oct2dec(13),m); mask = mask(1:end-1);
for ii=1:(m-1)
    % trellis state in binary
    bstate = de2bi(final_state,m-1);
    % xor'd feedback path in binary, fed into encoder input
    in = mod(sum(mask.*bstate),2);
    [tail,final_state] = convenc(in,trellis,final_state);
    x = [x; tail'];
end

end