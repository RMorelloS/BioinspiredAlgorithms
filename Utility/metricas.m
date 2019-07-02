%calculo das métricas utilizadas (dice, jaccard e dkl)
function [dice,jaccard,dklr] = metricas(Lims, sliceNII, imGroundTruth)
  imSegmentada = psrMultiLimiarizacao2(sliceNII, Lims, 1);
  LimsGT = unique(imGroundTruth);
  jaccard = feval('JAC', imSegmentada, LimsGT, imGroundTruth);
  dice = 2*jaccard/(1+jaccard);
  dklr = DKL(imSegmentada, imGroundTruth);
end
