function y = ldpc_enc(x,H,A,B,T,C,D,E)
n = size(H,2); k = n-size(H,1);
x = reshape(x,k,[]);
y = nan(n,size(x,2));
PHI = -E*inv(T)*B + D;
for ii=1:size(x,2)
    s = x(:,ii)';
    p1 = mod(-inv(PHI)*(-E*inv(T)*A + C)*s',2);
    p2 = mod(-inv(T)*(A*s' + B*p1),2);
    y(:,ii) = [s';p1;p2];
end
y = reshape(y,[],1);
end