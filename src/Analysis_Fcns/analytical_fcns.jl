
#############################
### Analytical functions for an Ising-X field in the zero hopping limit (t=0)
#############################
h4_analytic(Uu, beta, mum, L)=(4*Uu^2*(exp(beta*(Uu - mum))*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*
    (Uu^2 + mum^2)^(3/2)*(sqrt(Uu^2 + mum^2) + exp(2*beta*mum)*sqrt(Uu^2 + mum^2) + 4*exp(beta*(3*Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(beta*(Uu + mum - 2*sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) - 
    (2*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2) + 2*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + 2*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*
    (Uu^2 + mum^2) + exp(beta*(Uu + mum))*(Uu^2 + mum^2 - Uu*sqrt(Uu^2 + mum^2)) + exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(mum^2 + Uu*(Uu + sqrt(Uu^2 + mum^2))))^
    2/exp(2*beta*(-Uu + mum + 2*sqrt(Uu^2 + mum^2))) + (L^2*(2*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2) + 2*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + 
    2*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + exp(beta*(Uu + mum))*(Uu^2 + mum^2 - Uu*sqrt(Uu^2 + mum^2)) + 
    exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(mum^2 + Uu*(Uu + sqrt(Uu^2 + mum^2))))^2)/exp(2*beta*(-Uu + mum + 2*sqrt(Uu^2 + mum^2)))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*L^2*(Uu^2 + mum^2)^2)
;

h4_onsite_analytic(Uu, beta, mum, L)=(4*exp(beta*(Uu - mum))*Uu^2*(sqrt(Uu^2 + mum^2) + exp(2*beta*mum)*sqrt(Uu^2 + mum^2) + 4*exp(beta*(3*Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(beta*(Uu + mum - 2*sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*L^2*sqrt(Uu^2 + mum^2))
;

E_pot_analytic(Uu, beta, mum, L)=(-2*exp(Uu*beta)*Uu*(2*exp(3*Uu*beta) + 2/exp(beta*mum) + 2*exp(beta*mum) + exp(beta*(Uu - 2*sqrt(Uu^2 + mum^2)))*
    (1 - Uu/sqrt(Uu^2 + mum^2)) + exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2)))*(1 + Uu/sqrt(Uu^2 + mum^2))))/
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))
;

Sbil_B1_XX_analytic(Uu, beta, mum, L)=(16*(2*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2) + 2*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + 
    2*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + exp(beta*(Uu + mum))*(Uu^2 + mum^2 - Uu*sqrt(Uu^2 + mum^2)) + 
    exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(mum^2 + Uu*(Uu + sqrt(Uu^2 + mum^2))))^2)/
    (exp(2*beta*(-Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*L^2^2*(Uu^2 + mum^2)^2)
;
#added *L^2^2 below as compared to Mathematica file
χbil_B1_XX_analytic(Uu, beta, mum, L)=(8*((-2*exp(4*Uu*beta)*mum^2*(Uu^2 + mum^2))/Uu + 16*exp(2*Uu*beta)*beta*(Uu^2 + mum^2)^2 - (4*exp(2*beta*(Uu - mum))*(Uu^2 + mum^2)^2)/mum + 
    (4*exp(2*beta*(Uu + mum))*(Uu^2 + mum^2)^2)/mum - (8*exp(5*Uu*beta - beta*mum)*(2*Uu - mum)*(Uu^2 + mum^2)^2)/(Uu*mum) + 
    (8*exp(beta*(5*Uu + mum))*(2*Uu + mum)*(Uu^2 + mum^2)^2)/(Uu*mum) - (2*exp(8*Uu*beta)*(4*Uu^2 - mum^2)*(Uu^2 + mum^2)^2)/(Uu*mum^2) + 
    (exp(4*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2)*(4*Uu^3 + 3*Uu*mum^2 + 4*Uu^2*sqrt(Uu^2 + mum^2) + mum^2*sqrt(Uu^2 + mum^2)))/
    mum^2 - 2*exp(6*Uu*beta - 2*beta*sqrt(Uu^2 + mum^2))*(-4*Uu^4*beta - 6*Uu^2*beta*mum^2 - 2*beta*mum^4 + 4*Uu^3*beta*sqrt(Uu^2 + mum^2) + 
    mum^2*sqrt(Uu^2 + mum^2) + 4*Uu*beta*mum^2*sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(3*Uu + sqrt(Uu^2 + mum^2)))*
    (4*Uu^4*beta + 6*Uu^2*beta*mum^2 + 2*beta*mum^4 + 4*Uu^3*beta*sqrt(Uu^2 + mum^2) + mum^2*sqrt(Uu^2 + mum^2) + 
    4*Uu*beta*mum^2*sqrt(Uu^2 + mum^2)) + (8*exp(beta*(3*Uu + mum - 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2)*
    (-mum^2 + Uu*(-Uu + sqrt(Uu^2 + mum^2))))/(Uu - mum + sqrt(Uu^2 + mum^2)) + 
    (8*(Uu^2 + mum^2)*(-mum^2 + Uu*(-Uu + sqrt(Uu^2 + mum^2))))/(exp(beta*(-3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*
    (Uu + mum + sqrt(Uu^2 + mum^2))) + (exp(4*beta*(Uu - sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2)*
    (-mum^2 + 2*Uu*(-Uu + sqrt(Uu^2 + mum^2))))/(Uu + sqrt(Uu^2 + mum^2)) + 
    (8*exp(beta*(3*Uu - mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2)*(mum^2 + Uu*(Uu + sqrt(Uu^2 + mum^2))))/
    (-Uu - mum + sqrt(Uu^2 + mum^2)) + (8*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2)*
    (mum^2 + Uu*(Uu + sqrt(Uu^2 + mum^2))))/(-Uu + mum + sqrt(Uu^2 + mum^2))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(L^2)^2*(Uu^2 + mum^2)^2)
;

Sbil2_B1_XX_analytic(Uu, beta, mum, L)=(256*(-36*L^2*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^4 + 
    9*(L^2)^2*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^4 + 
    12*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*L^2*sqrt(Uu^2 + mum^2)*
    (2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2*
    (4*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) + 
    4*(-3*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2 + 
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2)*
    (4*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))))^2))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^4*(L^2)^6*(Uu^2 + mum^2)^2)
;

UL_bil_B1_analytic(Uu, beta, mum, L)=1-Sbil2_B1_XX_analytic(Uu, beta, mum, L)/(3*Sbil_B1_XX_analytic(Uu, beta, mum, L)^2);

Sbil_A1P_XX_analytic(Uu, beta, mum, L)=(16*(2*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2) + 2*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + 
    2*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + exp(beta*(Uu + mum))*
    (Uu^2 + mum^2 - Uu*sqrt(Uu^2 + mum^2)) + exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*
    (mum^2 + Uu*(Uu + sqrt(Uu^2 + mum^2))))^2)/(exp(2*beta*(-Uu + mum + 2*sqrt(Uu^2 + mum^2)))*
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*L^2^2*(Uu^2 + mum^2)^2) + 
    (16*(exp(beta*(Uu - mum))*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*(Uu^2 + mum^2)^(3/2)*(sqrt(Uu^2 + mum^2) + exp(2*beta*mum)*sqrt(Uu^2 + mum^2) + 
    4*exp(beta*(3*Uu + mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum - 2*sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 
    2*exp(beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) - 
    (3*(2*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2) + 2*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + 
    2*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + exp(beta*(Uu + mum))*
    (Uu^2 + mum^2 - Uu*sqrt(Uu^2 + mum^2)) + exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*
    (mum^2 + Uu*(Uu + sqrt(Uu^2 + mum^2))))^2)/exp(2*beta*(-Uu + mum + 2*sqrt(Uu^2 + mum^2)))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(L^2)^3*(Uu^2 + mum^2)^2) 
;

Sbil2_A1P_XX_analytic(Uu, beta, mum, L)=(64*(-2520*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^4 + 
    36*(L^2)^3*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^4 + 
    1680*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2)*
    (2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2*
    (4*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) - 
    140*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(Uu^2 + mum^2)*
    (4*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2 - 
    56*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(Uu^2 + mum^2)*
    (2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))*
    (16*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    8*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 8*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) + 
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^3*(Uu^2 + mum^2)^(3/2)*
    (64*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    32*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 32*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) - 
    168*(L^2)^2*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2*
    (3*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2 - 
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2)*
    (4*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))) + 
    4*L^2*(531*(2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^4 - 
    294*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2)*
    (2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2*
    (4*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    2*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) + 
    19*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(Uu^2 + mum^2)*
    (4*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 2*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*
    (-Uu + sqrt(Uu^2 + mum^2)) + 2*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))^2 + 
    6*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(Uu^2 + mum^2)*
    (2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2)))*
    (16*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2) + exp(beta*(Uu - mum))*sqrt(Uu^2 + mum^2) + exp(beta*(Uu + mum))*sqrt(Uu^2 + mum^2) + 
    8*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 8*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^4*(L^2)^7*(Uu^2 + mum^2)^2)
;

UL_bil_A1P_analytic(Uu, beta, mum, L)=1-Sbil2_A1P_XX_analytic(Uu, beta, mum, L)/(3*Sbil_A1P_XX_analytic(Uu, beta, mum, L)^2);


Sspin_A1_XX_analytic(Uu, beta, mum, L)= (4*exp(Uu*beta)*(2*exp(3*Uu*beta) + 2/exp(beta*mum) + 2*exp(beta*mum) + exp(beta*(Uu - 2*sqrt(Uu^2 + mum^2)))*(1 - Uu/sqrt(Uu^2 + mum^2)) + 
    exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2)))*(1 + Uu/sqrt(Uu^2 + mum^2))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*L^2)
;
#added below * (2/L^2) compared to Mathematica file
χspin_A1_XX_analytic(Uu, beta, mum, L)= 2*(2*exp(beta*(Uu - mum))*(-4*exp(beta*(3*Uu + mum))*Uu*(Uu^2 + mum^2) - 2*mum*(Uu^2 + mum^2) + 2*exp(2*beta*mum)*mum*(Uu^2 + mum^2) + 
    exp(beta*(Uu + mum - 2*sqrt(Uu^2 + mum^2)))*(2*Uu^3 + 2*Uu*mum^2 - 2*Uu^2*sqrt(Uu^2 + mum^2) - mum^2*sqrt(Uu^2 + mum^2)) + 
    exp(beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(2*Uu^3 + 2*Uu*mum^2 + 2*Uu^2*sqrt(Uu^2 + mum^2) + mum^2*sqrt(Uu^2 + mum^2))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*mum^2*(Uu^2 + mum^2) *L^2)
;

Sspin2_A1_XX_analytic(Uu, beta, mum, L)=(32*(2*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2) + 2*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + 
    2*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + exp(beta*(Uu + mum))*(Uu^2 + mum^2 - Uu*sqrt(Uu^2 + mum^2)) + 
    exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2 + Uu*sqrt(Uu^2 + mum^2)))^2)/
    (exp(2*beta*(-Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*L^2^2*(Uu^2 + mum^2)^2) + 
    (16*(exp(beta*(Uu - mum))*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*(Uu^2 + mum^2)^(3/2)*(sqrt(Uu^2 + mum^2) + exp(2*beta*mum)*sqrt(Uu^2 + mum^2) + 
    4*exp(beta*(3*Uu + mum))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + mum - 2*sqrt(Uu^2 + mum^2)))*(-Uu + sqrt(Uu^2 + mum^2)) + 
    2*exp(beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu + sqrt(Uu^2 + mum^2))) - 
    (3*(2*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2) + 2*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + 
    2*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2) + exp(beta*(Uu + mum))*(Uu^2 + mum^2 - Uu*sqrt(Uu^2 + mum^2)) + 
    exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2 + Uu*sqrt(Uu^2 + mum^2)))^2)/
    exp(2*beta*(-Uu + mum + 2*sqrt(Uu^2 + mum^2)))))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*L^2^3*(Uu^2 + mum^2)^2)
;

UL_magnetic_analytic(Uu, beta, mum, L)=1-Sspin2_A1_XX_analytic(Uu, beta, mum, L)/(3*Sspin_A1_XX_analytic(Uu, beta, mum, L)^2);


χ_charge_analytic(Uu, beta, mum, L)=(4*(2 + exp(beta*(Uu - mum)) + exp(beta*(Uu + mum)))*beta)/(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + 
    exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))
;
χ_pair_ZZ_analytic(Uu, beta, mum, L)=(8*exp(Uu*beta)*((2*sinh(beta*mum))/mum + (exp(Uu*beta)*sinh(2*beta*sqrt(Uu^2 + mum^2)))/sqrt(Uu^2 + mum^2)))/
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))
;
χ_pair_00_analytic(Uu, beta, mum, L)=(-4*(-4*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)))*Uu*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2)))*mum*
    sqrt(Uu^2 + mum^2) - 2*exp(beta*(Uu + 2*(mum + sqrt(Uu^2 + mum^2))))*mum*sqrt(Uu^2 + mum^2) + 
    exp(beta*(2*Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(-mum^2 + 2*Uu*(-Uu + sqrt(Uu^2 + mum^2))) + 
    exp(beta*(2*Uu + mum))*(mum^2 + 2*Uu*(Uu + sqrt(Uu^2 + mum^2)))))/(exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)))*
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*mum^2*sqrt(Uu^2 + mum^2))
;
heat_capacity_analytic(Uu, beta, mum, L)=(4*exp(beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(5*exp(2*beta*sqrt(Uu^2 + mum^2))*(Uu - mum)^2 + 
    16*exp(beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*mum^2 + exp(2*beta*(2*Uu + mum + sqrt(Uu^2 + mum^2)))*(-3*Uu + mum)^2 + 
    5*exp(2*beta*(mum + sqrt(Uu^2 + mum^2)))*(Uu + mum)^2 + exp(2*beta*(2*Uu + sqrt(Uu^2 + mum^2)))*(3*Uu + mum)^2 + 
    4*exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(6*Uu^2 + mum^2) + 
    5*exp(beta*(Uu + mum))*(2*Uu^2 + mum^2 - 2*Uu*sqrt(Uu^2 + mum^2)) + exp(beta*(5*Uu + mum + 4*sqrt(Uu^2 + mum^2)))*
    (2*Uu^2 + mum^2 - 2*Uu*sqrt(Uu^2 + mum^2)) + exp(beta*(5*Uu + mum))*(2*Uu^2 + mum^2 + 2*Uu*sqrt(Uu^2 + mum^2)) + 
    5*exp(beta*(Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(2*Uu^2 + mum^2 + 2*Uu*sqrt(Uu^2 + mum^2)) + 
    exp(2*Uu*beta)*(5*Uu^2 + 2*Uu*mum + 5*mum^2 - 4*Uu*sqrt(Uu^2 + mum^2) - 4*mum*sqrt(Uu^2 + mum^2)) + 
    exp(2*beta*(Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(5*Uu^2 - 2*Uu*mum + 5*mum^2 + 4*Uu*sqrt(Uu^2 + mum^2) - 
    4*mum*sqrt(Uu^2 + mum^2)) + exp(2*beta*(Uu + mum))*(5*Uu^2 - 2*Uu*(mum + 2*sqrt(Uu^2 + mum^2)) + 
    mum*(5*mum + 4*sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + 2*sqrt(Uu^2 + mum^2)))*
    (5*Uu^2 + 2*Uu*(mum + 2*sqrt(Uu^2 + mum^2)) + mum*(5*mum + 4*sqrt(Uu^2 + mum^2)))))/
    ((exp(beta*(2*Uu + mum)) + 4*exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2))) + 5*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2))) + 
    exp(beta*(4*Uu + mum + 2*sqrt(Uu^2 + mum^2))) + exp(beta*(2*Uu + mum + 4*sqrt(Uu^2 + mum^2))) + 
    4*exp(beta*(Uu + 2*(mum + sqrt(Uu^2 + mum^2)))))^2*L^2)
;

χ_charge_XX_analytic(Uu, beta, mum, L)=(2*(4*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)))*Uu*sqrt(Uu^2 + mum^2) - 2*exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2)))*mum*sqrt(Uu^2 + mum^2) + 
    2*exp(beta*(Uu + 2*(mum + sqrt(Uu^2 + mum^2))))*mum*sqrt(Uu^2 + mum^2) + exp(beta*(2*Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(2*Uu^2 + mum^2 - 2*Uu*sqrt(Uu^2 + mum^2)) - 
    exp(beta*(2*Uu + mum))*(2*Uu^2 + mum^2 + 2*Uu*sqrt(Uu^2 + mum^2))))/(exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)))*
    (5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*(L^2)*beta*mum^2*
    sqrt(Uu^2 + mum^2))
;


χ_charge_B1_analytic(Uu, beta, mum, L)=((-25*(Uu^2 + mum^2)^(3/2))/Uu + (exp(8*Uu*beta)*(Uu^2 + mum^2)^(3/2))/(3*Uu) + (2*exp(4*Uu*beta)*sqrt(Uu^2 + mum^2)*(6*Uu^2 + 7*mum^2))/Uu - 
    (4*exp(2*beta*(Uu - mum))*sqrt(Uu^2 + mum^2)*(Uu^4 - 14*Uu^3*mum - 19*Uu^2*mum^2 - 14*Uu*mum^3 - 20*mum^4))/(Uu^2*(3*Uu + 4*mum)) - 
    (10*exp(2*beta*(Uu - sqrt(Uu^2 + mum^2)))*(Uu^4 + 3*Uu^2*mum^2 + 2*mum^4))/Uu^2 + (10*exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))*(Uu^4 + 3*Uu^2*mum^2 + 2*mum^4))/Uu^2 + 
    (8*exp(2*Uu*beta)*sqrt(Uu^2 + mum^2)*(15*Uu^4 + 31*Uu^2*mum^2 + 16*mum^4))/(9*Uu^3 + 12*Uu*mum^2) + 
    (exp(5*Uu*beta - beta*mum)*(Uu^2 + mum^2)^(3/2)*(9*mum + 2*Uu*(-1 + beta*mum)))/(2*Uu*mum) + (exp(beta*(5*Uu + mum))*(Uu^2 + mum^2)^(3/2)*(9*mum + 2*Uu*(1 + beta*mum)))/
    (2*Uu*mum) + (5*exp(beta*(Uu - mum))*(Uu^2 + mum^2)^(3/2)*(-5*mum + 2*Uu*(-1 + 5*beta*mum)))/(2*Uu*mum) + 
    (5*exp(beta*(Uu + mum))*(Uu^2 + mum^2)^(3/2)*(-5*mum + 2*Uu*(1 + 5*beta*mum)))/(2*Uu*mum) + (exp(4*beta*(Uu - sqrt(Uu^2 + mum^2)))*Uu^2*sqrt(Uu^2 + mum^2))/
    (Uu - 2*sqrt(Uu^2 + mum^2)) + (exp(4*beta*(Uu + sqrt(Uu^2 + mum^2)))*Uu^2*sqrt(Uu^2 + mum^2))/(Uu + 2*sqrt(Uu^2 + mum^2)) - 
    (2*exp(6*Uu*beta - 2*beta*sqrt(Uu^2 + mum^2))*(Uu^2 + mum^2)*(Uu^2 + 2*mum^2 - 2*Uu*sqrt(Uu^2 + mum^2)))/(Uu*(5*Uu - 4*sqrt(Uu^2 + mum^2))) + 
    (2*exp(2*beta*(3*Uu + sqrt(Uu^2 + mum^2)))*(Uu^2 + mum^2)*(Uu^2 + 2*mum^2 + 2*Uu*sqrt(Uu^2 + mum^2)))/(Uu*(5*Uu + 4*sqrt(Uu^2 + mum^2))) - 
    (4*exp(2*beta*(Uu + mum))*(Uu^4 + 14*Uu^3*mum - 19*Uu^2*mum^2 + 14*Uu*mum^3 - 20*mum^4)*(Uu^2 + mum*(mum + sqrt(Uu^2 + mum^2))))/
    (Uu^2*(3*Uu - 4*mum)*(mum + sqrt(Uu^2 + mum^2))) - (exp(beta*(3*Uu - mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^4*(-3 + beta*mum) + 4*beta*mum^4*(mum - sqrt(Uu^2 + mum^2)) + 
    2*Uu^3*(mum + sqrt(Uu^2 + mum^2)) + 2*Uu*mum^2*(mum + sqrt(Uu^2 + mum^2)) + Uu^2*mum^2*(-3 + 5*beta*mum - 3*beta*sqrt(Uu^2 + mum^2))))/
    (mum*(mum - sqrt(Uu^2 + mum^2))) - (exp(beta*(3*Uu + mum - 2*sqrt(Uu^2 + mum^2)))*(Uu^4*(3 + beta*mum) + 4*beta*mum^4*(mum - sqrt(Uu^2 + mum^2)) + 
    2*Uu^3*(mum + sqrt(Uu^2 + mum^2)) + 2*Uu*mum^2*(mum + sqrt(Uu^2 + mum^2)) + Uu^2*mum^2*(3 + 5*beta*mum - 3*beta*sqrt(Uu^2 + mum^2))))/
    (mum*(mum - sqrt(Uu^2 + mum^2))) + (Uu^4*(-3 + beta*mum) + 2*Uu*mum^2*(mum - sqrt(Uu^2 + mum^2)) - 2*Uu^3*(-mum + sqrt(Uu^2 + mum^2)) + 
    4*beta*mum^4*(mum + sqrt(Uu^2 + mum^2)) + Uu^2*mum^2*(-3 + 5*beta*mum + 3*beta*sqrt(Uu^2 + mum^2)))/(exp(beta*(-3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*mum*
    (mum + sqrt(Uu^2 + mum^2))) + (exp(beta*(3*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*(Uu^4*(3 + beta*mum) + 2*Uu*mum^2*(mum - sqrt(Uu^2 + mum^2)) - 
    2*Uu^3*(-mum + sqrt(Uu^2 + mum^2)) + 4*beta*mum^4*(mum + sqrt(Uu^2 + mum^2)) + Uu^2*mum^2*(3 + 5*beta*mum + 3*beta*sqrt(Uu^2 + mum^2))))/
    (mum*(mum + sqrt(Uu^2 + mum^2))) + exp(3*Uu*beta - beta*mum)*(1 + exp(2*beta*mum))*Uu^2*sinh(2*beta*sqrt(Uu^2 + mum^2)))/
    ((5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(L^2)*beta*
    (Uu^2 + mum^2)^(3/2))
;



S_charge_B1_analytic(Uu, beta, mum, L)=((6*exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2)))*sqrt(Uu^2 + mum^2) + 5*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)))*sqrt(Uu^2 + mum^2) + 
    exp(beta*(4*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*sqrt(Uu^2 + mum^2) + 2*exp(beta*(Uu + 2*(mum + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2) + 
    exp(beta*(2*Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(-mum + sqrt(Uu^2 + mum^2)) + exp(beta*(2*Uu + mum))*(mum + sqrt(Uu^2 + mum^2)))*
    (2*exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2)))*sqrt(Uu^2 + mum^2) + 5*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)))*sqrt(Uu^2 + mum^2) + 
    exp(beta*(4*Uu + mum + 2*sqrt(Uu^2 + mum^2)))*sqrt(Uu^2 + mum^2) + 6*exp(beta*(Uu + 2*(mum + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2) + 
    exp(beta*(2*Uu + mum))*(-mum + sqrt(Uu^2 + mum^2)) + exp(beta*(2*Uu + mum + 4*sqrt(Uu^2 + mum^2)))*(mum + sqrt(Uu^2 + mum^2))))/
    (exp(2*beta*(mum + 2*sqrt(Uu^2 + mum^2)))*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))^2*(L^2)*(Uu^2 + mum^2))
;


χ_charge_A1P_analytic(Uu, beta, mum, L)=(4*exp(2*beta*sqrt(Uu^2 + mum^2))*(exp(Uu*beta) + 2*exp(beta*mum) + exp(beta*(Uu + 2*mum))))/
    ((exp(beta*(2*Uu + mum)) + 4*exp(beta*(Uu + 2*sqrt(Uu^2 + mum^2))) + 5*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2))) + 
    exp(beta*(4*Uu + mum + 2*sqrt(Uu^2 + mum^2))) + exp(beta*(2*Uu + mum + 4*sqrt(Uu^2 + mum^2))) + 
    4*exp(beta*(Uu + 2*(mum + sqrt(Uu^2 + mum^2)))))*L^2)
;



Gτ0_red_1_analytic(Uu, beta, mum, tau)=(exp(-(beta*(mum + 2*sqrt(Uu^2 + mum^2))) - (3*Uu + mum + 2*sqrt(Uu^2 + mum^2))*tau)*(exp(4*Uu*beta + beta*(mum + 2*sqrt(Uu^2 + mum^2)) + 2*(mum + sqrt(Uu^2 + mum^2))*tau)*
    sqrt(Uu^2 + mum^2) + 5*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)) + 2*(2*Uu + mum + sqrt(Uu^2 + mum^2))*tau)*sqrt(Uu^2 + mum^2) + 
    5*exp(2*beta*sqrt(Uu^2 + mum^2) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + Uu*(beta + 2*tau))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*sqrt(Uu^2 + mum^2) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + Uu*(beta + 6*tau))*sqrt(Uu^2 + mum^2) + exp(beta*(mum + 4*sqrt(Uu^2 + mum^2)) + 2*Uu*(beta + tau))*
    (-mum + sqrt(Uu^2 + mum^2)) + exp(2*beta*(mum + sqrt(Uu^2 + mum^2)) + Uu*(beta + 4*tau))*(-mum + sqrt(Uu^2 + mum^2)) + 
    exp(beta*mum + 4*sqrt(Uu^2 + mum^2)*tau + 2*Uu*(beta + tau))*(mum + sqrt(Uu^2 + mum^2)) + exp(2*beta*(mum + sqrt(Uu^2 + mum^2)) + 4*sqrt(Uu^2 + mum^2)*tau + Uu*(beta + 4*tau))*
    (mum + sqrt(Uu^2 + mum^2))))/(2*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + 
    exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2))
;
Gτ0_red_2_analytic(Uu, beta, mum, tau)=0
Gτ0_red_3_analytic(Uu, beta, mum, tau)=0
Gτ0_red_4_analytic(Uu, beta, mum, tau)=(exp(4*Uu*beta) + 5*exp(4*Uu*tau) + 5*exp(Uu*beta + beta*mum + 2*Uu*tau) + exp(Uu*beta + beta*mum + 6*Uu*tau) + 
    exp(-2*beta*sqrt(Uu^2 + mum^2) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + 2*Uu*(beta + tau))*(1 - mum/sqrt(Uu^2 + mum^2)) + 
    exp(-(beta*mum) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + Uu*(beta + 4*tau))*(1 - mum/sqrt(Uu^2 + mum^2)) + exp(Uu*beta - beta*mum + 4*Uu*tau + 2*mum*tau - 2*sqrt(Uu^2 + mum^2)*tau)*
    (1 + mum/sqrt(Uu^2 + mum^2)) + exp(2*(beta*sqrt(Uu^2 + mum^2) + mum*tau - sqrt(Uu^2 + mum^2)*tau + Uu*(beta + tau)))*(1 + mum/sqrt(Uu^2 + mum^2)))/
    (2*exp((3*Uu + mum)*tau)*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))))
;


G0τ_red_1_analytic(Uu, beta, mum, tau)=-(exp(4*Uu*beta) + 5*exp(4*Uu*tau) + 5*exp(Uu*beta + beta*mum + 2*Uu*tau) + exp(Uu*beta + beta*mum + 6*Uu*tau) + 
    exp(-2*beta*sqrt(Uu^2 + mum^2) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + 2*Uu*(beta + tau))*(1 - mum/sqrt(Uu^2 + mum^2)) + 
    exp(-(beta*mum) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + Uu*(beta + 4*tau))*(1 - mum/sqrt(Uu^2 + mum^2)) + exp(Uu*beta - beta*mum + 4*Uu*tau + 2*mum*tau - 2*sqrt(Uu^2 + mum^2)*tau)*
    (1 + mum/sqrt(Uu^2 + mum^2)) + exp(2*(beta*sqrt(Uu^2 + mum^2) + mum*tau - sqrt(Uu^2 + mum^2)*tau + Uu*(beta + tau)))*(1 + mum/sqrt(Uu^2 + mum^2)))/
    (2*exp((3*Uu + mum)*tau)*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2)))))
;
G0τ_red_2_analytic(Uu, beta, mum, tau)=0;
G0τ_red_3_analytic(Uu, beta, mum, tau)=0;
G0τ_red_4_analytic(Uu, beta, mum, tau)=-(exp(-(beta*(mum + 2*sqrt(Uu^2 + mum^2))) - (3*Uu + mum + 2*sqrt(Uu^2 + mum^2))*tau)*
    (exp(4*Uu*beta + beta*(mum + 2*sqrt(Uu^2 + mum^2)) + 2*(mum + sqrt(Uu^2 + mum^2))*tau)*sqrt(Uu^2 + mum^2) + 
    5*exp(beta*(mum + 2*sqrt(Uu^2 + mum^2)) + 2*(2*Uu + mum + sqrt(Uu^2 + mum^2))*tau)*sqrt(Uu^2 + mum^2) + 
    5*exp(2*beta*sqrt(Uu^2 + mum^2) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + Uu*(beta + 2*tau))*sqrt(Uu^2 + mum^2) + 
    exp(2*beta*sqrt(Uu^2 + mum^2) + 2*(mum + sqrt(Uu^2 + mum^2))*tau + Uu*(beta + 6*tau))*sqrt(Uu^2 + mum^2) + exp(beta*(mum + 4*sqrt(Uu^2 + mum^2)) + 2*Uu*(beta + tau))*
    (-mum + sqrt(Uu^2 + mum^2)) + exp(2*beta*(mum + sqrt(Uu^2 + mum^2)) + Uu*(beta + 4*tau))*(-mum + sqrt(Uu^2 + mum^2)) + 
    exp(beta*mum + 4*sqrt(Uu^2 + mum^2)*tau + 2*Uu*(beta + tau))*(mum + sqrt(Uu^2 + mum^2)) + 
    exp(2*beta*(mum + sqrt(Uu^2 + mum^2)) + 4*sqrt(Uu^2 + mum^2)*tau + Uu*(beta + 4*tau))*(mum + sqrt(Uu^2 + mum^2))))/
    (2*(5 + exp(4*Uu*beta) + 4*exp(beta*(Uu - mum)) + 4*exp(beta*(Uu + mum)) + exp(2*beta*(Uu - sqrt(Uu^2 + mum^2))) + exp(2*beta*(Uu + sqrt(Uu^2 + mum^2))))*sqrt(Uu^2 + mum^2))
;



