
function f = ackley(X)

% Ackley function.
n = length(X);
a = 20; b = 0.2; c = 2*pi;
s1 = 0; s2 = 0;

for i=1:n;
   s1 = s1+X(i)^2;
   s2 = s2+cos(c*X(i));
end
f = -a*exp(-b*sqrt(1/n*s1))-exp(1/n*s2)+a+exp(1);

