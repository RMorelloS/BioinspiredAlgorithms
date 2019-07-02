%normalização das imagens para o range 0-255
function slice = normalize(sliceNII)

minI = min(min(sliceNII));
maxI = max(max(sliceNII));
slice = ((sliceNII - minI)* 255)/(maxI  - minI);

end
