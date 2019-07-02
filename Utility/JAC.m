
%calculo do indice de jaccard (não da distância de jaccard)
%distancia de jaccard = 1 - indice de jaccard
function mediaIdx = JAC(imSegmentada, LimsGT, imGT)
    
    %maior valor de indice para a força-bruta de cada limiar do
    %ground-truth
    maiorJaccard = 0;
    %limiares da imagem segmentada
    Lims = unique(imSegmentada);
    %vetor com os maiores indices de força-bruta para cada limiar do GT
    vetorMaioresIdx = zeros(size(Lims,1));
    %para cada limiar do GT e cada limiar da imagem segmentada
    for i=1:size(LimsGT, 1)
        im1 = (imGT == LimsGT(i));
        for j=1:size(Lims,1)
            %coletando a camada correspondente aquele limiar
            im2 = (imSegmentada == Lims(j));
            %interseção das imagens
            intersectImg = im1 & im2;
            %união das imagens
            unionImg = im1 | im2;
            numerador = sum(intersectImg(:));
            denominador = sum(unionImg(:));
            %indice de jaccard para o limiar j do slice e limiar i do GT
            jaccardIdx = numerador/denominador;
            %se for maior que o já armazenado para os limiares j
            if jaccardIdx > maiorJaccard
                maiorJaccard = jaccardIdx;
            end
        end
        %vetor com os maiores indices de jaccard
        vetorMaioresIdx(i) = maiorJaccard;
        maiorJaccard = 0; 
    end
    %soma para tirar a media
    somaIdx=0;
    for i=1:size(vetorMaioresIdx,2)
        somaIdx=somaIdx+vetorMaioresIdx(i);
    end
    %indice de jaccard medio para todos os limiares
    mediaIdx = somaIdx/size(vetorMaioresIdx,2);
 

end