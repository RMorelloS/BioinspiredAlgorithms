function s = psrAvaliacaoTsallis2inicial(histograma, q, lim)   
    i = 1;
    if size(histograma, 1) == 1
       histograma = histograma';
    end
    lim = [1 lim 256];
    s = AvaliacaoTsallisv2(histograma,lim, i, q);
   
end