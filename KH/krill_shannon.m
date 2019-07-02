

function f = krill_shannon(X,H)

if ~isreal(X(1))
    X
    sdsdsd
end

L = [0 X(1) X(2) 255];

t2 = 0; SE = 0; EM = 0;
for t=2:length(L);    
    t1 = round(t2) + 1;
    t2 = round(L(t)) + 1;
    H1 = H(t1:t2);
    T = sum(H1);
    if (t2 > t1)
      EM = EM + log(t2-t1+1);
    end
    if T
      S = 0;
      H1 = H1/sum(H1);
      for k=1:length(H1)
         if H1(k)
            S = S - H1(k) * log(H1(k));
         end
      end
    end
    SE = SE + S;
end

if isreal(SE)
  f = SE;
else
  f = 10000;  
end

%f = (length(L) - 1)*EM - SE;

