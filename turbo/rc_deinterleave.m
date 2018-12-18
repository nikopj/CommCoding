function x = rc_deinterleave(x,R,C)
x = reshape(x',C,R)';
x = reshape(x,1,[]);
end