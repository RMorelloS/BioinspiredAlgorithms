
function ClanCenter =  CaculateClanCenter( OPTIONS, Clan, cindex)

ClanCenter = zeros(1, OPTIONS.numVar);

for Elephantindex = 1 : OPTIONS.numElephantInEachClan
    ClanCenter = ClanCenter + Clan{cindex}(Elephantindex).chrom;
end

ClanCenter = (1/OPTIONS.numElephantInEachClan)* ClanCenter;
