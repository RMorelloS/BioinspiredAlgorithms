function s = recursive_systems_evaluation(histograma,lim, i, q)
    s1 = TsallisEntropy(histograma(round(lim(i)):round(lim(i+1)-1)), q);
    if( i == size(lim, 2) - 2 )    
        s2 = TsallisEntropy(histograma(round(lim(i+1)):size(histograma,1)), q );
    else
        s2 = recursive_systems_evaluation(histograma,lim, i+1, q);   
    end 
    
    if( isinf(s1) || isinf(s2) )
        s = 0;
        return
    end
    s = s1+s2 + (1-q)*s1*s2;
end
