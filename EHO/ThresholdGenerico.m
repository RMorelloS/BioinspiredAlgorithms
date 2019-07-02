function [g,T1,T2,T3, T4] = ThresholdGenerico(L,e)
close all
% media e desvio-padrao da primeira gaussiana
m1 = (1/5)*L;
d1 = 7;
% media e desvio-padrao da segunda gaussiana
m2 = (2/5)*L;
d2 = 7;
% media e desvio-padrao da segunda gaussiana
m3 = (3/5)*L;
d3 = 7;
% media e desvio-padrao da segunda gaussiana
m4 = (4/5)*L;
d4 = 7;
% essas sao as duas gaussianas
g1 =  exp(-(((1:e:L)-m1).^2)/(2*d1^2))/(sqrt(2*pi)*d1);
g2 =  exp(-(((1:e:L)-m2).^2)/(2*d2^2))/(sqrt(2*pi)*d2);
g3 =  exp(-(((1:e:L)-m3).^2)/(2*d3^2))/(sqrt(2*pi)*d3);
g4 =  exp(-(((1:e:L)-m4).^2)/(2*d4^2))/(sqrt(2*pi)*d4);

% toma-se o valor maximo das ambas as gaussianas para normalizar adiante
[v1,T1] = max(g1);
[v2,T2] = max(g2);
[v3,T3] = max(g3);
[v4,T4] = max(g4);

% constroi-se uma unica gaussiana normaizada pelo maior valor
g = (g1 + g2 + g3+ g4)/sum(g1+g2+g3+g4);

