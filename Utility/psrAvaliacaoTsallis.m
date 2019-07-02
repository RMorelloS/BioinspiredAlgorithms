
function light=psrAvaliacaoTsallis(histograma, q, elemento)
n = size(elemento,2) + 2;
elemento=[1 elemento 256];
a = elemento(1);
b = elemento(2);
light = TsallisEntropy(histograma,q,a,b);
if isnan(light) light = 0; end
Plight = light;
for i=2:n-1
  a = elemento(i)+1;
  b = elemento(i+1);
  ES = TsallisEntropy(histograma,q,a,b);
  if ~isnan(ES)
      light = light + ES;
      Plight = Plight * ES;
  end
end

light = light + (1-q) * Plight;
end


function S = TsallisEntropy(histograma,q,a,b)
  a = round(a);
  b = round(b);
 if b > a
      H = histograma(a:b-1);
      H = H/sum(H);
      L = size(H,2);
      ret = sum(H.^q);
      S = (1-ret)/(q-1);
 else
     S = 0;
 end
end







