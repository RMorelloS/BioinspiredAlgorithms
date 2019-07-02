%função de preenchimento da struct
function resultado = preencher_struct(parametros,...
                                      i,q,j,dice,jaccard,dklr,cost,nomeArquivo,...
                                      t,thresholds,funcao,geracoes, Lims, idIndividuo, numGeracoes, melhores_avaliacoes)
    resultado.parametros = parametros;
    resultado.slice = i;
    resultado.q = q;
    resultado.iteracaoAlg = j;
    resultado.algoritmo = funcao;
    resultado.num_geracoes = geracoes;
    resultado.num_limiares = thresholds;
    resultado.dice = dice;
    resultado.jaccard = jaccard;
    resultado.dkl = dklr;
    resultado.entropia = cost;
    resultado.nome_arquivo = nomeArquivo;
    resultado.tempo = t;
    resultado.limiares = Lims;
    resultado.id_paciente = idIndividuo;
    resultado.numGeracoes = numGeracoes;                                    
    resultado.melhores_avaliacoes = melhores_avaliacoes;                                    
end
