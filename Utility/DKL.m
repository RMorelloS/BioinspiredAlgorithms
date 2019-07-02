%função de métrica: dkl
function DKLR = DKL(imSegmentada, imGroundTruth)
    hist1 = psrGrayHistogram(imSegmentada)';
    hist2 = psrGrayHistogram(imGroundTruth)';  
    distancia1 = 0;
    distancia2 = 0;
    for i=1:length(hist1)
        if hist1(i) && ~isnan(log(hist2(i)/hist1(i))) && ~isinf(log(hist2(i)/hist1(i)))
           distancia2 = distancia2 + (hist2(i) * log(hist2(i)/hist1(i)));
        end
        if hist2(i) && ~isnan(log(hist1(i)/hist2(i))) && ~isinf(log(hist1(i)/hist2(i)))
           distancia1 = distancia1 + (hist1(i) * log(hist1(i)/hist2(i)));
        end
    end
    DKLR = distancia1 + distancia2;
end

