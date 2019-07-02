%função de leitura das demais imagens
function V = leituraNII(path, nomeArquivo)
%caminho do arquivo 
pathNII = char(strcat(path, (strcat('\', nomeArquivo))));
%cabeçalho das imagens nii
info = nii_read_header(pathNII);
%matriz contendo todas as imagens do arquivo nii
V = nii_read_volume(info);
end
