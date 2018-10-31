function [G, H] = myBCHgen(m,t)
% inputs: m (order of galois field), t (error correcting capability)
% function returns the Generator and Parity check matricies of a BCH code
% defined by m and t. Matlab default primitive polynomials are used to 
% generate GF(2^m)

% m and t must satisfy these conditions
if (m<3 || t>2^(m-1))
    disp('error: invalid inputs');
    return;
end
n = 2^m-1;

% a primitive element in GF(2^m)
a = gf(2,m);
jj = 1:2:(2*t-1);
% vector of minimum polynomials
phi_vec = sym(zeros(1,length(jj)));
i=1;
for j=jj % population of phi_vec
    mp = minpol(a^j);
    phi_vec(i) = poly2sym( cast(mp.x,'like',1) );
    i=i+1;
end
% obtain generator polynomial through least common multiple
% and some acrobatics with the symbolic toolbox
lcm_sym = lcm(phi_vec);
lcm_poly = mod(sym2poly(lcm_sym), 2);
g = gf(lcm_poly,m);
k = n - length(g) + 1;

if k<1 || (m*t < n-k)
    fprintf('\nERROR\n');
    return
end

% making generator matrix, G, from circular shifts of generator polynomial
% vector, g
pad = zeros(1, mod(n-length(g), n) );
gx = [fliplr(g.x) pad];
G = zeros(k,n)*nan;
for i=1:k
    gx = circshift(gx,1);
    G(i,:) = gx;
end

% parity check matrix
H = zeros(length(jj)*m, n)*nan;
i=1;
for j=jj
    b = 0:(n-1);
    H(i:(i+m-1),:) = gftuple((b.*j)',m)';
    i=i+m;
end

% puts G in systematic form, G = [I P], where I is the kxk identity matrix
% making decoding much easier
G = mod(rref(G),2);

end