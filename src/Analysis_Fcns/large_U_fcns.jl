using LinearAlgebra


function Z0(beta, Uu, muM)
    return 4 + 2*exp(4*beta*Uu) + 4*exp(2*beta*(-muM + Uu)) + 
    4*exp(2*beta*(muM + Uu)) + exp(4*beta*(Uu - sqrt(muM^2 + Uu^2))) + 
    exp(4*beta*(Uu + sqrt(muM^2 + Uu^2)))
end

function χ_charge(beta, Uu, muM)
    return (4*beta*(2 + exp(2*beta*(-muM + Uu)) + 
        exp(2*beta*(muM + Uu))))/Z0(beta, Uu, muM)
end
function S_charge(beta, Uu, muM)
    return χ_charge(beta, Uu, muM)/beta
end

function χ_spin_XX(beta, Uu, muM)
    return ((-2*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*muM*Uu*sqrt(muM^2 + Uu^2) - 
    2*exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2)^(3/2) + 2*exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^(3/2) - 
    exp(2*beta*(muM + Uu))*muM*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*muM*
     (muM^2 + Uu*(Uu + sqrt(muM^2 + Uu^2)))))/(exp(2*beta*(muM - Uu + 2*sqrt(muM^2 + Uu^2)))*muM*(muM^2 + Uu^2)^(3/2)*Z0(beta, Uu, muM))
end


function S_spin_XX(beta, Uu, muM)
    return (2*exp(2*beta*Uu)*(2/exp(2*beta*muM) + 2*exp(2*beta*muM) + 2*exp(2*beta*Uu) + exp(2*beta*(Uu - 2*sqrt(muM^2 + Uu^2)))*
    (1 - Uu/sqrt(muM^2 + Uu^2)) + exp(2*beta*(Uu + 2*sqrt(muM^2 + Uu^2)))*(1 + Uu/sqrt(muM^2 + Uu^2))))/Z0(beta, Uu, muM)
end

function S_pair_ZZ(beta, Uu, muM)
    return (4*(2 + 2*exp(2*beta*(-muM + Uu)) + 2*exp(2*beta*(muM + Uu)) + exp(4*beta*(Uu - sqrt(muM^2 + Uu^2)))*(1 - Uu/sqrt(muM^2 + Uu^2)) + 
    exp(4*beta*(Uu + sqrt(muM^2 + Uu^2)))*(1 + Uu/sqrt(muM^2 + Uu^2))))/Z0(beta, Uu, muM)
end
function S_pair_00(beta, Uu, muM)
    return (4*(2 + 2*exp(2*beta*(-muM + Uu)) + 2*exp(2*beta*(muM + Uu)) + exp(4*beta*(Uu + sqrt(muM^2 + Uu^2)))*(1 - Uu/sqrt(muM^2 + Uu^2)) + 
    exp(4*beta*(Uu - sqrt(muM^2 + Uu^2)))*(1 + Uu/sqrt(muM^2 + Uu^2))))/Z0(beta, Uu, muM)
end
function χ_pair_ZZ(beta, Uu, muM)
    return (-2*exp(2*beta*Uu)*(2/exp(2*beta*muM) - 2*exp(2*beta*muM) + (exp(2*beta*(Uu - 2*sqrt(muM^2 + Uu^2)))*muM)/sqrt(muM^2 + Uu^2) - 
    (exp(2*beta*(Uu + 2*sqrt(muM^2 + Uu^2)))*muM)/sqrt(muM^2 + Uu^2)))/(muM*Z0(beta, Uu, muM))
end
function χ_pair_00(beta, Uu, muM)
    return (-4*exp(2*beta*(Uu + 2*sqrt(muM^2 + Uu^2)))*muM*sqrt(muM^2 + Uu^2) + 4*exp(2*beta*(Uu + 2*(muM + sqrt(muM^2 + Uu^2))))*muM*
    sqrt(muM^2 + Uu^2) + 8*exp(2*beta*(muM + 2*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2) + 
   exp(2*beta*(muM + 2*Uu + 4*sqrt(muM^2 + Uu^2)))*(2*muM^2 + 4*Uu*(Uu - sqrt(muM^2 + Uu^2))) - 
   2*exp(2*beta*(muM + 2*Uu))*(muM^2 + 2*Uu*(Uu + sqrt(muM^2 + Uu^2))))/(exp(2*beta*(muM + 2*sqrt(muM^2 + Uu^2)))*muM^2*
   sqrt(muM^2 + Uu^2)*Z0(beta, Uu, muM))  
end


function S_B1_nem_XY(beta, Uu, muM, Ls)
    return 1/Ls^4 *(32*(2*exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2) + 2*exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 
    2*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + exp(2*beta*(muM + Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + 
    exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu*(Uu + sqrt(muM^2 + Uu^2))))^2)/(exp(4*beta*(muM - Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^2*Z0(beta, Uu, muM)^2)
end
function S_A1p_Bil_XY(beta, Uu, muM, Ls)
    return S_B1_nem_XY(beta, Uu, muM, Ls) +
    1/Ls^6 *(32*(exp(2*beta*(muM - Uu + 2*sqrt(muM^2 + Uu^2)))*(4 + 2*exp(4*beta*Uu) + 4*exp(2*beta*(-muM + Uu)) + 4*exp(2*beta*(muM + Uu)) + exp(4*beta*(Uu - sqrt(muM^2 + Uu^2))) + 
    exp(4*beta*(Uu + sqrt(muM^2 + Uu^2))))*(muM^2 + Uu^2)^(3/2)*(exp(4*beta*sqrt(muM^2 + Uu^2))*sqrt(muM^2 + Uu^2) + 
    exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*sqrt(muM^2 + Uu^2) + 2*exp(2*beta*(muM + Uu))*(-Uu + sqrt(muM^2 + Uu^2)) + 
    2*exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*(Uu + sqrt(muM^2 + Uu^2))) + exp(2*beta*(muM - Uu + 2*sqrt(muM^2 + Uu^2)))*
   (4 + 2*exp(4*beta*Uu) + 4*exp(2*beta*(-muM + Uu)) + 4*exp(2*beta*(muM + Uu)) + exp(4*beta*(Uu - sqrt(muM^2 + Uu^2))) + exp(4*beta*(Uu + sqrt(muM^2 + Uu^2))))*
   (muM^2 + Uu^2)^(3/2)*(exp(4*beta*sqrt(muM^2 + Uu^2))*sqrt(muM^2 + Uu^2) + exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*sqrt(muM^2 + Uu^2) + 
    4*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*sqrt(muM^2 + Uu^2) + 2*exp(2*beta*(muM + Uu))*(-Uu + sqrt(muM^2 + Uu^2)) + 
    2*exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*(Uu + sqrt(muM^2 + Uu^2))) - 
  4*(2*exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2) + 2*exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 2*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*
      (muM^2 + Uu^2) + exp(2*beta*(muM + Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*
      (muM^2 + Uu*(Uu + sqrt(muM^2 + Uu^2))))^2))/(exp(4*beta*(muM - Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^2*Z0(beta, Uu, muM)^2)

end
function S_B1p_Bil_full(beta, Uu, muM, Ls)
    return S_B1_nem_XY(beta, Uu, muM, Ls) +
    1/Ls^4 *(64*exp(2*beta*(Uu - 2*(muM + 2*sqrt(muM^2 + Uu^2))))*(2*exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2) + 2*exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 
    2*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + exp(2*beta*(muM + Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + 
    exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu*(Uu + sqrt(muM^2 + Uu^2))))*(2*exp(2*beta*(muM + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 
    2*exp(2*beta*(Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 2*exp(2*beta*(Uu + 2*(muM + sqrt(muM^2 + Uu^2))))*(muM^2 + Uu^2) + 
    exp(2*beta*(muM + 2*Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + exp(2*beta*(muM + 2*Uu + 4*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu*(Uu + sqrt(muM^2 + Uu^2)))))/
  ((muM^2 + Uu^2)^2*Z0(beta, Uu, muM)^2)

end

function S_B1_nem_Z(beta, Uu, muM, Ls)
    return 1/Ls^4 *(16*(2*exp(2*beta*(muM + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 2*exp(2*beta*(Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 
    2*exp(2*beta*(Uu + 2*(muM + sqrt(muM^2 + Uu^2))))*(muM^2 + Uu^2) + exp(2*beta*(muM + 2*Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + 
    exp(2*beta*(muM + 2*Uu + 4*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu*(Uu + sqrt(muM^2 + Uu^2))))^2)/(exp(4*beta*(muM + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^2*Z0(beta, Uu, muM)^2)
end

function χ_B1_nem_XY(beta, Uu, muM, Ls)
    return 1/Ls^4 *(32*(-(exp(8*beta*Uu)*(-4*beta*muM^2 + Uu)*(muM^2 + Uu^2)) + 8*beta*exp(4*beta*Uu)*(muM^2 + Uu^2)^2 - (exp(4*beta*(-muM + Uu))*(muM^2 + Uu^2)^2)/muM + 
    (exp(4*beta*(muM + Uu))*(muM^2 + Uu^2)^2)/muM - (2*exp(-2*beta*muM + 6*beta*Uu)*(-muM + 2*Uu)*(muM^2 + Uu^2)^2)/Uu^2 + 
    (exp(2*beta*(muM + 3*Uu))*(muM^2 + Uu^2)^2*(-muM - 4*Uu + sqrt(muM^2 + Uu^2)))/Uu^2 - 
    (exp(8*beta*(Uu - sqrt(muM^2 + Uu^2)))*(-2*muM^2*Uu - 2*Uu^3 + muM^2*sqrt(muM^2 + Uu^2) + 2*Uu^2*sqrt(muM^2 + Uu^2)))/4 + 
    (exp(8*beta*(Uu + sqrt(muM^2 + Uu^2)))*(2*muM^2*Uu + 2*Uu^3 + muM^2*sqrt(muM^2 + Uu^2) + 2*Uu^2*sqrt(muM^2 + Uu^2)))/4 + 
    exp(8*beta*Uu - 4*beta*sqrt(muM^2 + Uu^2))*(2*beta*muM^4 + 6*beta*muM^2*Uu^2 + 4*beta*Uu^4 - (muM^2*sqrt(muM^2 + Uu^2))/2 - 4*beta*muM^2*Uu*sqrt(muM^2 + Uu^2) - 
      4*beta*Uu^3*sqrt(muM^2 + Uu^2)) + exp(4*beta*(2*Uu + sqrt(muM^2 + Uu^2)))*(2*beta*muM^4 + 6*beta*muM^2*Uu^2 + 4*beta*Uu^4 + (muM^2*sqrt(muM^2 + Uu^2))/2 + 
      4*beta*muM^2*Uu*sqrt(muM^2 + Uu^2) + 4*beta*Uu^3*sqrt(muM^2 + Uu^2)) + 
    ((muM^2 + Uu^2)*(2*Uu^3 + muM^2*(muM - sqrt(muM^2 + Uu^2)) + 2*muM*Uu*(muM - sqrt(muM^2 + Uu^2)) + Uu^2*(muM - sqrt(muM^2 + Uu^2))))/
     (exp(2*beta*(muM - 3*Uu + 2*sqrt(muM^2 + Uu^2)))*Uu^2) + (exp(2*beta*(muM + 3*Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*
      (2*Uu^3 + 2*muM*Uu*(muM - sqrt(muM^2 + Uu^2)) + muM^2*(-muM + sqrt(muM^2 + Uu^2)) + Uu^2*(-muM + sqrt(muM^2 + Uu^2))))/Uu^2 + 
    (2*exp(2*beta*(muM + 3*Uu - 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*(Uu^2 + muM*(muM + sqrt(muM^2 + Uu^2))))/Uu + 
    (exp(-2*beta*muM + 6*beta*Uu + 4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2)*(2*Uu^3 + muM^2*(muM + sqrt(muM^2 + Uu^2)) + 2*muM*Uu*(muM + sqrt(muM^2 + Uu^2)) + 
       Uu^2*(muM + sqrt(muM^2 + Uu^2))))/Uu^2 - (2*(muM^2 + Uu^2)^2*cosh(2*beta*sqrt(muM^2 + Uu^2))*(exp(2*beta*muM)*(muM + sqrt(muM^2 + Uu^2)) - 
       2*(-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*sqrt(muM^2 + Uu^2)*cosh(2*beta*muM) + 2*(1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM*sinh(2*beta*muM)))/
     (exp(2*beta*(-3*Uu + sqrt(muM^2 + Uu^2)))*Uu^2)))/((muM^2 + Uu^2)^2*Z0(beta, Uu, muM)^2)
end
function χ_A1p_Bil_XY(beta, Uu, muM, Ls)
    return χ_B1_nem_XY(beta, Uu, muM, Ls)+
    1/Ls^6 *(8*((-4*beta*(2*exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2) + 2*exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 
    2*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + exp(2*beta*(muM + Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + 
    exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))^2)/(exp(4*beta*(muM - Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^2) + 
((4 + 2*exp(4*beta*Uu) + 4*exp(2*beta*(-muM + Uu)) + 4*exp(2*beta*(muM + Uu)) + exp(4*beta*(Uu - sqrt(muM^2 + Uu^2))) + exp(4*beta*(Uu + sqrt(muM^2 + Uu^2))))*
  (4*beta*exp(2*beta*(Uu + 6*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^(3/2) + 4*beta*exp(2*beta*(2*muM + Uu + 6*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^(3/2) + 
   exp(2*beta*(muM + 2*Uu + 4*sqrt(muM^2 + Uu^2)))*(-8*beta*muM^2*Uu - 8*beta*Uu^3 + 8*beta*Uu^2*sqrt(muM^2 + Uu^2) + muM^2*(-1 + 4*beta*sqrt(muM^2 + Uu^2))) + 
   exp(2*beta*(muM + 2*Uu + 8*sqrt(muM^2 + Uu^2)))*(8*beta*muM^2*Uu + 8*beta*Uu^3 + 8*beta*Uu^2*sqrt(muM^2 + Uu^2) + muM^2*(1 + 4*beta*sqrt(muM^2 + Uu^2)))))/
 (exp(2*beta*(muM + 6*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)^(3/2)) + 
(4*(-((beta*(2*exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2) + 2*exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + 
        2*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2) + exp(2*beta*(muM + Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + 
        exp(2*beta*(muM + Uu + 4*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))^2)/exp(4*beta*(muM - Uu + 2*sqrt(muM^2 + Uu^2)))) - 
   2*(8*beta*exp(4*beta*Uu)*muM^4 + 4*beta*exp(8*beta*Uu)*muM^4 + 2*beta*exp(4*beta*(2*Uu + sqrt(muM^2 + Uu^2)))*muM^4 + 
     2*beta*exp(8*beta*Uu - 4*beta*sqrt(muM^2 + Uu^2))*muM^4 + 16*beta*exp(4*beta*Uu)*muM^2*Uu^2 + 4*beta*exp(8*beta*Uu)*muM^2*Uu^2 + 
     6*beta*exp(4*beta*(2*Uu + sqrt(muM^2 + Uu^2)))*muM^2*Uu^2 + 6*beta*exp(8*beta*Uu - 4*beta*sqrt(muM^2 + Uu^2))*muM^2*Uu^2 + 8*beta*exp(4*beta*Uu)*Uu^4 + 
     4*beta*exp(4*beta*(2*Uu + sqrt(muM^2 + Uu^2)))*Uu^4 + 4*beta*exp(8*beta*Uu - 4*beta*sqrt(muM^2 + Uu^2))*Uu^4 + 4*beta*exp(4*beta*(2*Uu + sqrt(muM^2 + Uu^2)))*
      muM^2*Uu*sqrt(muM^2 + Uu^2) - 4*beta*exp(8*beta*Uu - 4*beta*sqrt(muM^2 + Uu^2))*muM^2*Uu*sqrt(muM^2 + Uu^2) + 4*beta*exp(4*beta*(2*Uu + sqrt(muM^2 + Uu^2)))*
      Uu^3*sqrt(muM^2 + Uu^2) - 4*beta*exp(8*beta*Uu - 4*beta*sqrt(muM^2 + Uu^2))*Uu^3*sqrt(muM^2 + Uu^2) - (exp(4*beta*(-muM + Uu))*(muM^2 + Uu^2)^2)/muM + 
     (exp(4*beta*(muM + Uu))*(muM^2 + Uu^2)^2)/muM + (exp(-2*beta*muM + 6*beta*Uu)*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 
        (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(-muM + sqrt(muM^2 + Uu^2)) - 
     (exp(2*beta*(muM + 3*Uu - 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
        (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(-muM + sqrt(muM^2 + Uu^2)) + 
     (exp(2*beta*(muM + 3*Uu))*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
        (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(muM + sqrt(muM^2 + Uu^2)) - 
     ((muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*
         sqrt(muM^2 + Uu^2)))/(exp(2*beta*(muM - 3*Uu + 2*sqrt(muM^2 + Uu^2)))*(muM + sqrt(muM^2 + Uu^2))) + 
     (exp(-2*beta*muM + 6*beta*Uu)*(muM^2 + Uu^2)*(muM - exp(4*beta*muM)*muM + (1 + exp(4*beta*muM))*sqrt(muM^2 + Uu^2))*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 
        (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/Uu^2 - 
     (exp(6*beta*Uu - 2*beta*(muM + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*((-1 + exp(4*beta*muM))*muM + (1 + exp(4*beta*muM))*sqrt(muM^2 + Uu^2))*
       ((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/Uu^2 + 
     (exp(8*beta*Uu)*sqrt(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))^2*muM^2 + 2*(1 + exp(8*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
        2*(-1 + exp(8*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/4 - (exp(8*beta*(Uu - sqrt(muM^2 + Uu^2)))*sqrt(muM^2 + Uu^2)*
       ((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))^2*muM^2 + 2*(1 + exp(8*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 2*(-1 + exp(8*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/
      4) + (4 + 2*exp(4*beta*Uu) + 4*exp(2*beta*(-muM + Uu)) + 4*exp(2*beta*(muM + Uu)) + exp(4*beta*(Uu - sqrt(muM^2 + Uu^2))) + exp(4*beta*(Uu + sqrt(muM^2 + Uu^2))))*
    sqrt(muM^2 + Uu^2)*(4*beta*exp(4*beta*Uu)*(muM^2 + Uu^2)^(3/2) + beta*exp(2*beta*(-muM + Uu))*(muM^2 + Uu^2)^(3/2) + 
     beta*exp(2*beta*(muM + Uu))*(muM^2 + Uu^2)^(3/2) + exp(4*beta*(Uu - sqrt(muM^2 + Uu^2)))*(-2*beta*muM^2*Uu - 2*beta*Uu^3 + 2*beta*Uu^2*sqrt(muM^2 + Uu^2) + 
       (muM^2*(-1 + 4*beta*sqrt(muM^2 + Uu^2)))/4) + (exp(4*beta*(Uu + sqrt(muM^2 + Uu^2)))*(8*beta*muM^2*Uu + 8*beta*Uu^3 + 8*beta*Uu^2*sqrt(muM^2 + Uu^2) + 
        muM^2*(1 + 4*beta*sqrt(muM^2 + Uu^2))))/4)))/(muM^2 + Uu^2)^2))/Z0(beta, Uu, muM)^2
end

function χ_B1p_Bil_full(beta, Uu, muM, Ls)
    return χ_B1_nem_XY(beta, Uu, muM, Ls)+
    1/Ls^4 *(16*(32*beta*exp(4*beta*Uu)*muM^4 + 64*beta*exp(4*beta*Uu)*muM^2*Uu^2 + 32*beta*exp(4*beta*Uu)*Uu^4 - (2*exp(4*beta*(-muM + Uu))*(muM^2 + Uu^2)^2)/muM + 
    (2*exp(4*beta*(muM + Uu))*(muM^2 + Uu^2)^2)/muM - (exp(4*beta*(Uu - sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 2*Uu^2 - 
       2*Uu*sqrt(muM^2 + Uu^2)))/Uu + (exp(8*beta*Uu - 4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 2*Uu^2 - 
       2*Uu*sqrt(muM^2 + Uu^2)))/Uu + (exp(4*beta*Uu)*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 2*Uu^2 - 2*Uu*sqrt(muM^2 + Uu^2)))/
     (-Uu + 2*sqrt(muM^2 + Uu^2)) - (exp(8*beta*(Uu - sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 2*Uu^2 - 
       2*Uu*sqrt(muM^2 + Uu^2)))/(-Uu + 2*sqrt(muM^2 + Uu^2)) - (2*exp(2*beta*(-muM + Uu))*(muM^2 + Uu^2)*(-muM^2 - Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/
     (-muM - Uu + sqrt(muM^2 + Uu^2)) + (2*exp(2*beta*(muM + 3*Uu - 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*(-muM^2 - Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/
     (-muM - Uu + sqrt(muM^2 + Uu^2)) - (2*exp(2*beta*(muM + Uu))*(muM^2 + Uu^2)*(-muM^2 - Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(muM - Uu + sqrt(muM^2 + Uu^2)) + 
    (4*(muM^2 + Uu^2)*(-muM^2 - Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(exp(2*beta*(muM - 3*Uu + 2*sqrt(muM^2 + Uu^2)))*(muM - Uu + sqrt(muM^2 + Uu^2))) - 
    (2*exp(2*beta*(muM + Uu))*(muM^2 + Uu^2)*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(-muM + Uu + sqrt(muM^2 + Uu^2)) + 
    (2*exp(-2*beta*muM + 6*beta*Uu + 4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2)*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(-muM + Uu + sqrt(muM^2 + Uu^2)) - 
    (2*exp(2*beta*(-muM + Uu))*(muM^2 + Uu^2)*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(muM + Uu + sqrt(muM^2 + Uu^2)) + 
    (2*exp(2*beta*(muM + 3*Uu + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(muM + Uu + sqrt(muM^2 + Uu^2)) + 
    (2*exp(-2*beta*muM + 6*beta*Uu)*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
       (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(-muM + sqrt(muM^2 + Uu^2)) - 
    (2*exp(2*beta*(muM + 3*Uu - 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
       (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(-muM + sqrt(muM^2 + Uu^2)) + 
    (2*exp(2*beta*(muM + 3*Uu))*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
       (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(muM + sqrt(muM^2 + Uu^2)) - 
    (2*(muM^2 + Uu^2)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*
        sqrt(muM^2 + Uu^2)))/(exp(2*beta*(muM - 3*Uu + 2*sqrt(muM^2 + Uu^2)))*(muM + sqrt(muM^2 + Uu^2))) + 
    (2*exp(2*beta*(muM + 2*Uu - 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*(-(exp(2*beta*(muM + 2*sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*(-muM + sqrt(muM^2 + Uu^2))) + 
       exp(2*beta*Uu)*muM*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + exp(2*beta*(Uu + 2*sqrt(muM^2 + Uu^2)))*muM*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2))))/
     (muM*(muM - sqrt(muM^2 + Uu^2))) - (2*(muM^2 + Uu^2)*(exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2)*(muM + sqrt(muM^2 + Uu^2)) + 
       exp(2*beta*(muM + Uu))*muM*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) + exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*muM*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2))))/
     (exp(4*beta*(muM - Uu + sqrt(muM^2 + Uu^2)))*muM*(muM + sqrt(muM^2 + Uu^2))) + 
    ((muM^2 + Uu^2)*((exp(2*beta*(muM + 8*Uu))*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 2*Uu^2 - 2*Uu*sqrt(muM^2 + Uu^2)))/(-Uu + 2*sqrt(muM^2 + Uu^2)) + 
       (2*exp(14*beta*Uu)*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)))/(-muM - Uu + sqrt(muM^2 + Uu^2)) + 
       (2*exp(2*beta*(2*muM + 7*Uu))*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)))/(muM - Uu + sqrt(muM^2 + Uu^2)) + 
       (exp(2*beta*(muM + 8*Uu - 2*sqrt(muM^2 + Uu^2)))*(-((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2) - 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2)))/Uu))/
     exp(2*beta*(muM + 6*Uu)) + ((muM^2 + Uu^2)*((2*exp(2*beta*(2*muM + 7*Uu))*(-muM^2 - Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(-muM - Uu + sqrt(muM^2 + Uu^2)) + 
       (exp(2*beta*(muM + 8*Uu - 2*sqrt(muM^2 + Uu^2)))*(-((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2) - 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2)))/
        (-Uu + 2*sqrt(muM^2 + Uu^2))))/exp(2*beta*(muM + 4*Uu + 2*sqrt(muM^2 + Uu^2))) - 
    (exp(4*beta*Uu)*(muM^2 + Uu^2)*(muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/Uu + 
    (exp(8*beta*Uu)*(muM^2 + Uu^2)*(muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/Uu - 
    (exp(4*beta*(Uu - sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*(muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/
     (Uu + 2*sqrt(muM^2 + Uu^2)) + (exp(4*beta*(2*Uu + sqrt(muM^2 + Uu^2)))*(muM^2 + Uu^2)*
      (muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/(Uu + 2*sqrt(muM^2 + Uu^2)) + 
    ((muM^2 + Uu^2)*(((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + 2*Uu^2 - 2*Uu*sqrt(muM^2 + Uu^2))/(exp(4*beta*(-2*muM - 4*Uu + sqrt(muM^2 + Uu^2)))*Uu) + 
       (2*exp(6*beta*muM + 14*beta*Uu)*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
          (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(-muM + sqrt(muM^2 + Uu^2)) + 
       (2*exp(2*beta*(5*muM + 7*Uu))*((1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*muM^2 + (1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu^2 + 
          (-1 + exp(4*beta*sqrt(muM^2 + Uu^2)))*Uu*sqrt(muM^2 + Uu^2)))/(muM + sqrt(muM^2 + Uu^2)) - 
       (exp(8*beta*muM + 12*beta*Uu)*(muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/Uu))/exp(8*beta*(muM + Uu)) + 
    exp(-2*beta*muM - 8*beta*Uu + 4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + Uu^2)*((2*exp(14*beta*Uu)*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/
       (-muM + Uu + sqrt(muM^2 + Uu^2)) + (2*exp(2*beta*(2*muM + 7*Uu))*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(muM + Uu + sqrt(muM^2 + Uu^2)) + 
      (exp(2*beta*(muM + 8*Uu - 2*sqrt(muM^2 + Uu^2)))*(muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/Uu + 
      (exp(2*beta*(muM + 8*Uu))*(muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/(Uu + 2*sqrt(muM^2 + Uu^2))) + 
    ((muM^2 + Uu^2)*((-2*exp(2*beta*(2*muM + 7*Uu))*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(-muM + Uu + sqrt(muM^2 + Uu^2)) - 
       (2*exp(14*beta*Uu)*(muM^2 + Uu^2 + Uu*sqrt(muM^2 + Uu^2)))/(muM + Uu + sqrt(muM^2 + Uu^2)) - 
       (exp(2*beta*(muM + 8*Uu - 2*sqrt(muM^2 + Uu^2)))*(muM^2 + exp(4*beta*sqrt(muM^2 + Uu^2))*(muM^2 + 2*Uu^2 + 2*Uu*sqrt(muM^2 + Uu^2))))/
        (Uu + 2*sqrt(muM^2 + Uu^2))))/exp(2*beta*(muM + 6*Uu))))/((muM^2 + Uu^2)^2*Z0(beta, Uu, muM)^2)
end

function χ_B1_nem_Z(beta, Uu, muM, Ls)
    return 1/Ls^4 *(4*(2*exp(4*beta*(muM + 2*Uu + sqrt(muM^2 + Uu^2)))*muM^3*sqrt(muM^2 + Uu^2) - 2*exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*muM*(muM^2 + Uu^2)^(3/2) - 
    8*exp(2*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*muM*(muM^2 + Uu^2)^(3/2) - 8*exp(2*beta*(3*muM + Uu + 2*sqrt(muM^2 + Uu^2)))*muM*(muM^2 + Uu^2)^(3/2) - 
    4*exp(4*beta*(Uu + sqrt(muM^2 + Uu^2)))*Uu*(muM^2 + Uu^2)^(3/2) + 4*exp(4*beta*(2*muM + Uu + sqrt(muM^2 + Uu^2)))*Uu*(muM^2 + Uu^2)^(3/2) + 
    32*beta*exp(4*beta*(muM + Uu + sqrt(muM^2 + Uu^2)))*muM*Uu*(muM^2 + Uu^2)^(3/2) + 4*exp(2*beta*(muM + 3*Uu))*muM*(muM^2 + Uu^2)*(-muM - Uu + sqrt(muM^2 + Uu^2)) + 
    4*exp(6*beta*(muM + Uu))*muM*(muM^2 + Uu^2)*(muM - Uu + sqrt(muM^2 + Uu^2)) + 4*exp(6*beta*(muM + Uu) + 8*beta*sqrt(muM^2 + Uu^2))*muM*(muM^2 + Uu^2)*
     (-muM + Uu + sqrt(muM^2 + Uu^2)) + 4*exp(2*beta*(muM + 3*Uu + 4*sqrt(muM^2 + Uu^2)))*muM*(muM^2 + Uu^2)*(muM + Uu + sqrt(muM^2 + Uu^2)) - 
    exp(4*beta*(muM + 2*Uu - sqrt(muM^2 + Uu^2)))*muM*Uu*(muM^2 + Uu^2 - Uu*sqrt(muM^2 + Uu^2)) - 
    2*exp(4*beta*(muM + Uu))*muM*Uu*(muM^2 - 4*beta*muM^2*sqrt(muM^2 + Uu^2) - 8*beta*Uu^2*sqrt(muM^2 + Uu^2) + 8*beta*Uu*(muM^2 + Uu^2)) + 
    exp(4*beta*(muM + 2*Uu + 3*sqrt(muM^2 + Uu^2)))*muM*Uu*(muM^2 + Uu*(Uu + sqrt(muM^2 + Uu^2))) + 2*exp(4*beta*(muM + Uu + 2*sqrt(muM^2 + Uu^2)))*muM*Uu*
     (muM^2 + 4*beta*(2*muM^2*Uu + 2*Uu^3 + muM^2*sqrt(muM^2 + Uu^2) + 2*Uu^2*sqrt(muM^2 + Uu^2)))))/
  (exp(4*beta*(muM + sqrt(muM^2 + Uu^2)))*muM*Uu*(muM^2 + Uu^2)^(3/2)*Z0(beta, Uu, muM)^2)
end