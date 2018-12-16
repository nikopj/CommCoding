function y = ldpc_msgpass2(x,H,M,MCi,LCi,MVj,LVj)
n=size(H,2); k=n-size(H,1);
x = reshape(x,n,[]);


y = nan(size(x));
% mij, message node i, check node j
for kk=1:size(x,2)
    xx = x(:,kk);
    llr0 = log((1-M(2,xx))./(M(2,xx)))';
    llr0(llr0==Inf)=300; llr0(llr0==-Inf)=-300;
    llr = llr0;
    % initialize m matrix
    mij = repmat(llr,1,n-k);
    mji = mij';
    % defining max iterations
    inds = ones(1,n);
    for iters=1:1:3
       if sum(mod(H*(llr<0),2))==0
           break;
       end
       P = nan(1,n-k);
       for ii=1:n
           if inds(ii)==0
               continue;
           end
           Ci = MCi(ii,1:LCi(ii));
           S = sum(mji(Ci,:));
           for jj=1:(n-k)
               % ci -> fj
               mij(ii,jj) = llr0(ii) + S(ii) - mji(jj,ii);
               % fj -> ci
               Vj = MVj(jj,1:LVj(jj));
               if sum(Vj==ii)
                   P(jj)= prod(tanh(mij(Vj,jj)/2));
               end
               t = P(jj) / tanh(mij(ii,jj)/2);
               mji(jj,ii) = log((1+t)/(1-t));
           end
           llr(ii) = llr0(ii) + sum(mji(Ci,ii));
           if abs(llr(ii))>30
               inds(ii)=0;
           end
       end
       
    end
    y(:,kk) = llr<0;
end
y = y(1:k,:);
y = reshape(y,[],1);
end