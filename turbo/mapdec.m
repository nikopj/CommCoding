function Le12 = mapdec(ys,yp,trellis,Lc,Le21,N)
num_states = trellis.numStates;
num_in = trellis.numInputSymbols; 
num_out= trellis.numOutputSymbols;
% map decoding
% STATES, S,SP
% S: state+input = next state
% SP: state+input = prev. state
S = trellis.nextStates+1; S_minus=S(:,2); S_plus=S(:,1);
SP = zeros(size(S));
for ii=1:num_states
    for jj=1:num_in
        SP(S(ii,jj),jj) = ii;
    end
end
SP_minus=SP(:,2); SP_plus=S(:,1);

% alpha and beta are really alpha,beta tilda
alpha = zeros(N+1,num_states); alpha(1,1) = 1;
beta = zeros(N+1,num_states); beta(N+1,1) = 1;
gamma = nan(N+1,num_states,num_states);
gamma_e = nan(N+1,num_states,num_states);
Le12 = zeros(N,1);

% parity bit generation (prev. state, next state)
x_parity = zeros(num_states,num_states);
uk = zeros(num_states,num_states);
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

Le21=[0,Le21];
ys=[0,ys]; yp=[0,yp];
for k=2:(N+1)
    % GAMMA
   yks = ys(k); ykp = yp(k);
   for s=1:size(S,1)
       for sp=SP(s,1:num_in)
           disp(sp)
           xkp = x_parity(sp,s);
           u   = uk(sp,s);
           gamma_e(k,sp,s)   = exp(0.5*Lc*ykp*xkp);
           gamma(k,sp,s)   = exp(0.5*u*(Le21(k)+Lc*yks))*gamma_e(k,sp,s);
       end
   end
   % ALPHA
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
            sp=SP_minus(t);
            den = den + alpha(k-1,sp)*gamma_e(k,sp,s)*beta(k,s);
        end
    end
    Le12(k-1) = log(num/(den));
end
end