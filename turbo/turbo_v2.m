% Nikola Janjusevic
% Turbo 2
close all, clear all;

num_run = 1; % number of runs
num_sym = 1e3; % number of symbols
snr_vec = 20;

n = 2;
k = 1;
m = 5;
rate = k/n;

trellis = poly2trellis(m,[31,27],31);
msg = randi([0 1], num_sym, 1);
[x,final_state] = convenc(msg,trellis);
tail=x(end-1:end);

% forcing the RSCC to the zero state
% mask is the feedback path connections
mask = de2bi(oct2dec(31),m); mask = mask(1:end-1);
for ii=1:(m-1)
    % trellis state in binary
    bstate = de2bi(final_state,m-1);
    % xor'd feedback path in binary, fed into encoder input
    in = mod(sum(mask.*bstate),2);
    [tail,final_state] = convenc(in,trellis,final_state);
    x = [x; tail'];
end
% total of 2*(m-1) tailbits appended

tx = -2*x+1;
rx = tx + 10^(-snr_vec(1)/20)*randn(size(tx));

N = num_sym;
num_states = trellis.numStates;
num_in = trellis.numInputSymbols; 
num_out= trellis.numOutputSymbols;
% map decoding
% STATES, S,SP
% S: state+input = next state
% SP: state+input = prev. state
S = trellis.nextStates+1; S_minus=S(:,2); S_plus=S(:,1);
SP = nan(size(S));
for ii=1:num_states
    for jj=1:num_in
        SP(S(ii,jj),jj) = ii;
    end
end
SP_minus=SP(:,2); SP_plus=S(:,1);

% alpha and beta are really alpha,beta tilda
alpha = zeros(N+1,num_states); alpha(1,1) = 1;
beta = zeros(N+1,num_states); beta(N+1,1) = 1;
gamma = zeros(N+1,num_states,num_states);
gamma_e = zeros(N+1,num_states,num_states);
Le = zeros(N+1,1);
Lc = 2*rate*10^(snr_vec(1)/10);

% parity bit generation (prev. state, next state)
x_parity = nan(num_states,num_states);
uk = nan(num_states,num_states);
for ii=1:num_states
    for input=1:num_in
        enc = de2bi(trellis.outputs(ii,input),log2(num_out),'left-msb');
        par = bi2de(enc(2:end));
        jj  = S(ii,input); % next state
        x_parity(ii,jj) = -2*par+1;
        % input value to the encoder
        uk(ii,jj) = -2*(input-1)+1;
    end
end

y = reshape(rx,2,[]);
% GAMMA
for k=1:(N+1)
    yk = y(:,k); yks = yk(1); ykp = yk(2:end);
   for s=1:size(S,1)
       for sp=SP(s,1:num_in)
           xkp = x_parity(sp,s);
           u   = uk(sp,s);
           gamma_e(k,sp,s)   = exp(0.5*Lc*ykp*xkp);
           gamma(k,sp,s)   = exp(0.5*u*(Le(k)+Lc*yks))*gamma_e(k,sp,s);
       end
   end
end
% ALPHA
for k=2:(N+1)
   den = 0;
   for s=1:size(S,1)
       num=0;
       for sp=SP(s,1:num_in)
           num = num + alpha(k-1,sp)*gamma(k,sp,s);
           den = den + alpha(k-1,sp)*gamma(k,sp,s);
       end
       alpha(k,s) = num;
   end
   alpha(k,:) = alpha(k,:)/den;
end
% BETA
for k=(N+1):-1:3
    den = 0;
    for sp=1:size(SP,1)
       num=0;
       for s=S(sp,1:num_in)
           num = num + beta(k,s)*gamma(k,sp,s);
           den = den + alpha(k-1,sp)*gamma(k,sp,s);
       end
       beta(k-1,sp) = num;
    end
    beta(k-1,:) = beta(k-1,:)/den;
end
% Le12
for k=2:(N+1)
    num=0; den=0;
    for s=1:size(S,1)
        for t=1:size(SP_plus,1)
            sp=SP_plus(t);
            num = num + alpha(k-1,sp)*gamma_e(k,sp,s)*beta(k,s);
        end
        for t=1:size(SP_minus,1)
            sp=SP_plus(t);
            den = den + alpha(k-1,sp)*gamma_e(k,sp,s)*beta(k,s);
        end
    end
    Le(k-1) = Lc*y(1,k-1) + log(num/den);
end
msg_hat = Le(1:end-1)<0;
biterr(msg,msg_hat)

