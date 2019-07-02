%função para retornar parametros dos algoritmos bio-inspirados
function parametros = retorna_parametros(funcao)
    switch char(funcao)
        case 'KH'
            parametros.NK = 40;
            parametros.UB = 253;
            parametros.LB = 2;
        case 'EHO'
            parametros.popsize = 100;
            parametros.numClan = 5;
            parametros.Keep = 2;
            parametros.alpha = 0.5;
            parametros.beta = 0.1;
        case 'CS' 
            parametros.n = 40;
            parametros.pa = 0.5;
            parametros.UB = 253;
            parametros.LB = 2;
        case 'FF'
            parametros.nFireflies = 50;
            parametros.method = 'TE';

        case 'GWO'
            parametros.nLobos = 30;
            parametros.UB = 253;
            parametros.LB = 2;
    
        case 'GOA'
            parametros.nGafanhotos = 30;
            parametros.UB = 253;
            parametros.LB = 2;
    
        case 'WOA'
            parametros.UB = 253;
            parametros.LB = 2;
            parametros.nBaleias = 30;
    end
    
end
