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
# Enthropy [kJ/kmolK] TO DO

using Roots

# Number of components (1)Nitrogen - (2)Methane - (3)Ethane - (4)Propane - (5)Butane

const NC = 5 # Number of components

const Pc = [33.5, 45.4, 48.2, 41.9, 37.5 ]*1.e5  # [N/m2]
const Tc = [126.2, 190.6, 305.4, 369.8, 425.2 ]  # [K]
const w =  [0.040, 0.008, 0.098, 0.152, 0.193 ]  # acentric factor [-]

# Ideal has heat capacity (See Hid below):
const Cp = [ 7.440 -0.324e-2  6.400e-6 -2.790e-9; # Nitrogen
             4.598  1.245e-2  2.86e-6  -2.703e-9; # Methane
             1.292  4.254e-2 -1.657e-5  2.081e-9; # Ethane
            -1.009  7.315e-2 -3.789e-5  7.678e-9; # Propane
             2.266  7.913e-2 -2.647e-5 -0.674e-9] # n-butane

const Tref = 298.15  # for Ideal gas heat capacity

const kinteraction = zeros(NC,NC)     # Here: SRK binary interaction parameters set to zero

const R = 8.314 # [kJ/kmolK]

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
    Hid += x[i] * (t1 + t2 + t3 + t4)
  end
  Hid *= 4.184 #From kcal to kJ
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
