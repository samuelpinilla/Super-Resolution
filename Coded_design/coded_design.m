function [gdmd] = coded_design(Shots,N)
rep = 2;
for ind_corrida=1:8
    shots = Shots;
    M = N;
    tam_mask=2*shots+1;
    alargue=floor(tam_mask/2);
    gdmd = zeros(M+alargue*2,M+alargue*2,shots);
    for rrep=2:rep
        for sh=1:length(shots)
            gdmd_aux=gdmddesign(M,shots(sh));
            gdmd(alargue+1:alargue+M,alargue+1:alargue+M,:)=gdmd_aux;
            MH=0;
            w=[1 5 .5];
            for j=1:1
                gdmd3=changelocalvh(gdmd,alargue,alargue,MH,tam_mask);
            end
            gdmd=gdmd3(alargue+1:alargue+M,alargue+1:alargue+M,:);
        end
    end    
end
end