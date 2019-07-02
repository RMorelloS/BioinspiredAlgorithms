function vetor=ordenaVetor(vetor)

ultimo = vetor(end);
vetor = sort(vetor(1:end-1));

vetor = [vetor, ultimo];

end
