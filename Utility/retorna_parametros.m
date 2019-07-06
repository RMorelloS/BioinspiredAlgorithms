%função para retornar parametros dos algoritmos bio-inspirados
function parametros = retorna_parametros(funcao)
    switch char(funcao)
        case 'KH'
            parametros.pop_size = 40;
            parametros.UB = 253;
            parametros.LB = 2;
        case 'EHO'
            parametros.pop_size = 100;
            parametros.numClan = 5;
            parametros.Keep = 2;
            parametros.alpha = 0.5;
            parametros.beta = 0.1;
        case 'CS' 
            parametros.pop_size = 40;
            parametros.pa = 0.5;
            parametros.UB = 253;
            parametros.LB = 2;
        case 'FF'
            parametros.pop_size = 50;
            parameters.UB = 253;
            parameters.LB = 2;
        case 'GWO'
            parametros.pop_size = 30;
            parametros.UB = 253;
            parametros.LB = 2;
    
        case 'GOA'
            parametros.pop_size = 30;
            parametros.UB = 253;
            parametros.LB = 2;
    
        case 'WOA'
            parametros.UB = 253;
            parametros.LB = 2;
            parametros.pop_size = 30;
    end
    
end
