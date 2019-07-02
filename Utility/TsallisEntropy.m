
function S = TsallisEntropy(H,q)

      if(size(H,1) <= 1)
          S = - Inf;
          return
      end
      
      soma = sum(H);
      if soma == 0
          S = - Inf;
          return
      end
      H = H/soma;
      ret = sum(H.^q);
      
      S = (1-ret)/(q-1);
end
