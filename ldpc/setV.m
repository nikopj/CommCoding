% Vj set of all message nodes connected to check node fj
function S = setV(j,I,J)
    S = J(I==j);
end