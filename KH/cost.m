function f=cost(X,H, q)
%funcao custo de Ackley
%f = ackley(X);
%f = morefood(X);
%f = krill_shannon(X,H);
if q == 1
    f=feval('psrAvaliacaoShannon',H, X');
else
    f=feval('psrAvaliacaoTsallis2inicial',H, q, X');         
end

if isnan(f)
  f = 1000;
end

