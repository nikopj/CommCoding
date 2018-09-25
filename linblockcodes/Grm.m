function G = Grm(r,m)
if r==0
    G = ones(1,2^m);
    return
elseif r==m
    row = zeros(1,2^m); row(1,end) = 1;
    G = [ Grm(m-1,m) ; row ];
    return
else
    x = Grm(r-1,m-1);
    G = [ Grm(r,m-1) Grm(r,m-1) ; zeros(size(x)) x ];
    return
end
end