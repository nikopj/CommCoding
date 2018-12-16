% Ci set of all check nodes connected to message node ci
function S = setC(i,I,J)
    S = I(J==i);
end