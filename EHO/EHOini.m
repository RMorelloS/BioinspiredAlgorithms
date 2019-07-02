function EHOini
 
    geracoes = 100;
    thresholds = 3;
    q = 0.9;
    parametros.popsize = 200;
    parametros.numClan = 5;
    %quantos elefantes deixar para a próxima geração
    parametros.Keep = 2;
    %fator de influencia da matriarca no resto do clã
    parametros.alpha = 0.5;
    %fator de influencia do elefante com a melhor solução nos demais
    parametros.beta = 0.1;
    im1 = imread('C:\Users\Windows 10Pro\Desktop\Populars\teste (12).jpg');
    [Lims, entropia] = EHO(im1, thresholds, geracoes, q, parametros);
    Lims
    entropia

end