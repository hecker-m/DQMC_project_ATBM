Bravais_dirdir2dir(l::Lattice)=Bravais_dirdir2dir(Bravais(l))


function Bravais_dirdir2dir(B::Bravais)
    lat=B.l
    N=length(B)

    Bsrctrg2dir = lat[:Bravais_srctrg2dir]
    Bsrcdir2trg = lat[:Bravais_srcdir2trg]
    output = Array{Int}(undef, N, N)

    src1=1;
    for d1 in 1:N
        trg1=Bsrcdir2trg[src1, d1]
        for d2 in 1:N
            src2=trg1
            trg2=Bsrcdir2trg[src2, d2]
            d3=Bsrctrg2dir[src1, trg2]
            output[d1, d2]=d3
        end

    end
    return output

end



