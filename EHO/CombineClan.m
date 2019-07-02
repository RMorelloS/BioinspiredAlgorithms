
function Population = CombineClan(OPTIONS, NewClan)

j =1;
popindex = 1;
while popindex <= OPTIONS.popsize
    for clanindex = 1 : OPTIONS.numClan
        Population(popindex)   = NewClan{clanindex}(j);
        popindex = popindex + 1;
    end
    j = j+1;
end
