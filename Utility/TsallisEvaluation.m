function s = TsallisEvaluation(histograma, q, lim)   
    i = 1;
    if size(histograma, 1) == 1
       histograma = histograma';
    end
    lim = [1 lim 256];
    s = recursive_systems_evaluation(histograma,lim, i, q);
   
end