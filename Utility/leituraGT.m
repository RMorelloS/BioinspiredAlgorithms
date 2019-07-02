%função de leitura da imagem de ground-truth
function [V,individuo] = leituraGT(path)
%nome do individuo, começando sempre a partir da última \ 
%do caminho da pasta
split = strsplit(path,'\');
individuo = split(size(split, 2));
%caminho do arquivo da segmentação manual
arquivoGT = char((strcat('\', strcat(individuo, '_seg.nii'))));
groundtruth = char(strcat(path, arquivoGT));
%cabeçalho das imagens nii
info = nii_read_header(groundtruth);
%matriz contendo todas as imagens do arquivo nii
V = nii_read_volume(info);
end
