################################################################################
###                               SRK PACKAGE                                ###
### Based on Linhart & Skogestad, 2008                                       ###
### Corrected by Leguizamón, 2016                                            ###
### Translated to Julia by Leguizamón, 2017                                  ###
################################################################################

# The following document contains a list of functions to calculate different   #
# properties using the Soave-Redlich-Kwong equation of state.
# For a given Temperature[K], Pressure [Pa] and Composition, it is possible
# to obtain:
# Compressibility Factor (for each phase using Kamath's, 2010
#                         derivative criteria)
# Fugacity coefficients (for each phase)
# Enthalpy [kJ/kmol]
# Molar volume [m3/kmol]
# Enthropy [kJ/kmolK] 

using Roots

# Number of components (1)Nitrogen - (2)Methane - (3)Ethane - (4)Propane - (5)Butane

const NC = 5 # Number of components

const Pc = [33.5, 45.4, 48.2, 41.9, 37.5 ]*1.e5  # [N/m2]
const Tc = [126.2, 190.6, 305.4, 369.8, 425.2 ]  # [K]
const w =  [0.040, 0.008, 0.098, 0.152, 0.193 ]  # acentric factor [-]

const R = 8.314 # [kJ/kmolK]

# Ideal has heat capacity (See Hid below):
const Cp = [ 3.539 -0.261e-3 0.007e-5   0.157e-8 -0.099e-11;     # Nitrogen
             4.568 -8.975e-3 3.631e-5  -3.407e-8  1.091e-11;     # Methane
             4.178 -4.427e-3 5.660e-5  -6.651e-8  2.487e-11;     # Ethane
             3.847  5.131e-3 6.011e-5  -7.893e-8  3.079e-11;     # Propane
             5.547  5.536e-3 8.057e-5 -10.571e-8  4.134e-11] * R # n-butane

const Tref = 298.15  # for Ideal gas heat capacity

const kinteraction = zeros(NC,NC)     # Here: SRK binary interaction parameters set to zero



ZRA = 0.29056 - 0.08775 * w

## Validation of inputs
# Composition
function validateX(x)
  if sum(x) < 0 || sum(x) > 1 || any(x) < 0 || any(x) > 1 || len(x)!= NC
    error("Composition poorly specified")
  end
end
# Temperature
function validateT(T)
  if T <= 0
    error("Temperature badly specified")
  end
end
# Pressure
function validateP(P)
  if P <= 0
    error("Pressure badly specified")
  end
end

## Calculations for given T, P  and composition (x)

function getPar(T, P, x)
  Tre = T./Tc
  Pre = P./Pc
  m = 0.480 + 1.574.*w - 0.176.*w.^2 # Eq 15
  a = (1+m.*(1-Tre.^0.5)).^2; # Eq 13
  Ap = 0.42747.*a.*Pre./Tre.^2; # Eq 8
  Bp = 0.08664.*Pre./Tre # Eq 9

  Ab = zeros(NC, NC)
  for i = 1:NC
      for j = 1:NC
          Ab[i,j]=(Ap[i]*Ap[j])^0.5
      end
  end
  # Mixture a and b
  A = 0.0;
  for i = 1:NC
      for j = 1:NC
          A += x[i]*x[j]*Ab[i,j]*(1-kinteraction[i,j]);
      end
  end
  B = 0.0;
  for i=1:NC
      B += x[i]*Bp[i];
  end
  return A, B, Ap, Bp, m
end

## Compressibility factor

# Get all real solutions for the equation within 0 and 1
function getZs(A, B)
  f(y) = y^3 -y^2 +(A-B-B^2)*y -A*B
  z = fzeros(f, 0, 1)
end

# Return Z for a given phase
function phaseZ(z, phase)
  zOut = 0.0
  if phase == "liquid" || phase == 1
    zOut = min(z)
  elseif phase == "vapor" || phase == "vapour" || phase == 2
    zOut = max(z)
  else
    error("Bad specified phase")
  end
  return zOut
end

# Finds the phases in a given vector of compressibility factors
# using Kamath's derivative criteria
function findPhaseZ(z, A, B)
  firstDer(y) = 3.0*y^2 - 2.0*y + A-B-B^2
  secondDer(y) = 6.0 * y - 2.0
  zVap = 0.0
  zLiq = 0.0
  for i = 1:length(z)
     if firstDer(z[i]) >= 0 && secondDer(z[i]) >= 0
       zVap = z[i]
     elseif firstDer(z[i]) >= 0 && secondDer(z[i]) <= 0
       zLiq = z[i]
     end
   end
   return zVap, zLiq
 end


## Fugacity coefficient:
function getPhi(T, P, x, A, B, Ap, Bp, Z)
  # Calculates the fugacity coefficient for a given set of SRK parameters
  # and compressibility
  corrphi = zeros(NC)
  for i = 1:NC
     for j = 1:NC
       corrphi[i] += x[j] * Ap[j]^0.5 * (1.0-kinteraction[i,j])
     end
  end
  phi = zeros(NC)
  factor1 = 0.0
  for i in 1:NC
    factor1 = (Z-1.0)*Bp[i]/B-log(Z-B) - A/B*log((Z+B)/Z)
    phi[i] = exp(factor1*((corrphi[i]*2.0*Ap[i]^0.5/A)-Bp[i]/B)) # Eq 21 Soave,1972
  end
  return phi
end

## Enthalpy calculations
# Calculates ideal gas enthalpy
function getHideal(T, P, x)
  Hid = 0.0
  for i = 1:NC
    t1 = (Cp[i,1]) * (T - Tref)
    t2 = 1/2 * Cp[i,2] * (T^2 - Tref^2)
    t3 = 1/3 * Cp[i,3] * (T^3 - Tref^3)
    t4 = 1/4 * Cp[i,4] * (T^4 - Tref^4)
    t5 = 1/5 * Cp[i,5] * (T^5 - Tref^5)
    Hid += x[i] * (t1 + t2 + t3 + t4 + t5)
  end
  return Hid
end

# Departure enthalpy
function getdadT(T, P , x, m, Ap, A, B)
  # Calculates the factor dadT
  dadT = 0.0
  for i = 1:NC
    for j = 1:NC
      f1 = -R / 2.0 * sqrt(abs(0.42748 / T)) * x[i] * x[j]
      t1 = m[j] * sqrt(abs(Ap[i]*(Tc[j]/P*(T^2)/(Pc[j])*(R^2))))
      t2 = m[i] * sqrt(abs(Ap[j]*(Tc[i]/P*(T^2)/(Pc[i])*(R^2))))
      dadT -= f1 * (t1+t2)
    end
  end
  return dadT
end

# Calculates the departure enthalpy given dadT and Z
function getHSRK(T, P, x, A, B, dadT, Z)
  HSRK = 0.0
  f1 = 1/ (B*R*T/P)
  f2 = A*R^2*T^2/P - T*dadT
  f3 = log(Z / (Z + B))
  t1 = R * T * (Z - 1)
  HSRK = f1 * f2 * f3 + t1
  return HSRK
end

# Total enthalpy
function getHreal(Hideal, HSRK)
  # Wrapper function for the calculation of the real enthalpy
  Hreal = Hideal + HSRK
end

## Molar volume

# For Liquid Phase with Peneleoux correction

function getVLiq(T, P, x, Z)
  c=0
  for i=1:NC
    f1 = (R * Tc[i]) / (Pc[i])
    c += x[i] * (0.40768 * (0.29441 - ZRA[i]) * f1)
  end
V = ((Z * R * T / P)- c)
end
# For vapour phase
function getVVap(T, P, x, Z)
  V = Z * R * T / P
end

# Getting Molar volume
function getV(T, P, x, Z)
    if Z<0.1
        V = getVLiq(T, P, x, Z)
     else
        V = getVVap(T, P, x, Z)
    end
    return V
end

## Entropy

function getSSRK(T, P, x, Z)
    A, B, Ap, Bp, m = getPar(T, P, x)
    dadT = getdadT(T, P , x, m, Ap, A, B)
    corrb = (R*T/P)
    Vt = getV(T, P, x, Z)
    fact1 = dadT/ ( B * corrb )
    fact2 = log(Vt / ( Vt + B * corrb))
    Ssrk = fact1 * fact2 - R*log(Z * (1 - B * corrb/Vt))
    return Ssrk
end


# Ideal gas entropy

function getSideal(T, P, x)
  # Constant pressure contribution
  SidP = 0.0
  for i = 1:NC
    t1 = Cp[i,1] * log(T - Tref)
    t2 = Cp[i,2] * (T - Tref)
    t3 = 1/2 * Cp[i,3] * (T^2 - Tref^2)
    t4 = 1/3 * Cp[i,4] * (T^3 - Tref^3)
    t5 = 1/4 * Cp[i,5] * (T^4 - Tref^4)
    SidP = SidP + x[i] * (t1 + t2 + t3 + t4 + t5)
  end
  #   Constant temperature contribution
  SidT = R * log(P/Pref);

  #  Mixing contribution
  SidM = 0;
  for i =1:NC
      if x[i] > 0
          SidM = SidM + x[i] * log(x[i])
      end
  end
  S = SidP - SidT - SidM
  return SIdeal
end

# Final entropy
function getS(T, P, x, Z)
  SSrk = getSSRK(T, P, x, Z)
  SIdeal = getSideal(T, P, x)
  return SSRk + SIdeal
end

# Test
# Not updated for S
function TestSRK(T, P, x)
  A, B, Ap, Bp, m = getPar(T, P, x)
  comps = getZs(A, B)
  compV, compL = findPhaseZ(comps, A, B)
  phiV = getPhi(T, P, x, A, B, Ap, Bp, compV)
  phiL = getPhi(T, P, x, A, B, Ap, Bp, compL)
  Hid = getHideal(T, P, x)
  dadT = getdadT(T, P , x, m, Ap, A, B)
  HSRK = getHSRK(T, P, x, A, B, dadT, compV)
  Hreal = getHreal(Hid,HSRK)
  return compV, compL, phiV, phiL, Hid, Hreal
end

T = 303.15
P = 101325.0
x = [0.0, 0.0, 0.0, 0.0, 1.0]
@time print(TestSRK(T, P, x))
