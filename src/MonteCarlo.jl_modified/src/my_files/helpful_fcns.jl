"""
to_string(x::Real)
transforms a Float into a string that can be used for filenames.
    E.g. 10.034 → 10p034
"""
function to_string(x::Real) 
    x2=Int(floor(x*10^6))
    trigger=false;
    dig = 0;
    st="";
   while (x2 >0 && dig <7) || dig <7
       dig +=1;
       if dig ==7
           trigger==false ? st="p0" * st : st="p" * st;
           trigger=true
           #If trigger stayed zero till digit=7, there was no .digit, thus we write xxxp0,
           #otherwise xxxpyyy
       end
       m2=mod(x2,10)
       
       if trigger || m2!=0
            if dig ==7
                st=string(x2) * st
                #Once we are past the 7th digit, the whole integer number is added upfront
            else
                st=string(m2) * st
                trigger=true;   
                #trigger stays zero for as long as the right-most digits are zero 
                #Once it is true, all digits are recorded
            end
       end
       x2=div(x2, 10);
   end
   return st
end


function _zero!(A::AbstractArray{G}) where {G}
    for i in eachindex(A)
        if abs(A[i])<10^(-13)
            A[i]=zero(G)
        end
    end
    return A
end
