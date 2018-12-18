function x = rc_interleave(x,R,C)
x = reshape(x,R,C);
x = reshape(x',1,[]);
end
