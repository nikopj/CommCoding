function E = myBCHerror(S, ygf_enc, m, t)
% error vector generator following the 
% Peterson Direct-Solution Decoding algorithm
% inputs: 
%   ygf_enc: column wise BCH(m,t) encoded vector in Galois array
%   S: Syndrom matrix (v*Htranspose)
%   m: order of galois field
%   t: error correcting capability of code

E = zeros(size(ygf_enc));
alpha_vec = gf(2,m).^(1:(2*t-1));
% loop through the received vectors with errors 
for i=find(sum(S,1))
    y = ygf_enc(:,i);
    % received vector evaluated
    p = polyval(y,alpha_vec);
    px = p.x;
    A = zeros(t)*nan;
    k = 0;
    % A matrix construction
    for j=1:t
        pad_len = t - k;
        if pad_len > 0
            pad = [1 zeros(1,pad_len-1)];
        else
            pad = [];
        end
        A(j,:) = [ fliplr( px(max(1,k-t+1):k) ) pad ];
        k = k+2; % sydrom rows increase in length by two every iteration
    end
    A = gf(A,m);
    sz = t;
    L = [0];
    
    % Peterson method of solving A*Lambda = b
    while( det(A) == 0 && sz>1 )
        A = A(1:(end-2),1:(end-2));
        sz=sz-2;
    end
    if det(A) ~= 0
        b = p(1:2:(2*sz-1))';
        L = [A\b ; 1];
    end
    r = roots(L);
    
    % ugly checking for multiplicity of roots
    mult = false;
    for j=1:length(r)
        for k=1:length(r)
            if r(j) == r(k)
                mult = true;
                break;
            end
        end
    end
    
    if length(r) ~= t || mult
        % cannot decode, leave alone
        E(:,i) = zeros(size(E,1),1);
    else
        % error vector
        E(:,r.x) = 1;
    end
    
end
end