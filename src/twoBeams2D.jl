include("pf_SQMC.jl")
include("plotsAndTests.jl")
include("sobol.jl")

const Δt = 0.02
#The size of the random acceleration
const σa = 10.0 #m/s^2
#The size of the random change in ϕ
const σϕ = 0.5
#The size of the random change in β
const σβ = 0.005
#The size of the noise on acceleration measurements
const σam = .05
const σym = .01
#Wavelength in meters
const λ = 0.13

# State x = [pₓ, py, vₓ, vy, ϕ₁, ϕ₂, β₁, β₂, r₁, r₂, ax, ay]
# Noise u = [ax, ay, ϕ₁, ϕ₂, β₁, β₂, ny1, ny2, nax, nay]
# Out   y₁ = ∑rₖe^(iβₖ)e^(i2π<pₖ,u(ϕₖ)>)+ny
#       y₂ = a + na

const ixp, ixv, ixϕ, ixβ, ixr, ixa = 1:2, 3:4, 5:6, 7:8, 9:10, 11:12
const ina, inϕ, inβ, iny, inya = 1:2, 3:4, 5:6, 7:8, 9:10
const s1, s2 = [40, 40], [-40, 40]

#OBS LET NOISE BE NORMAL
@inbounds function signalUpdate!(x, t, n, xt, i)
  a = n[ina]
  xt[ixp] = x[ixp] + Δt*x[ixv] + Δt^2/2*a #p
  xt[ixv] = x[ixv] + Δt*a                 #v
  xt[ixϕ] = x[ixϕ] + Δt*σϕ*n[inϕ]         #ϕ
  xt[ixβ] = x[ixβ] - ((cos(xt[ixϕ])-(cos(x[ixϕ])))*x[ixp[1]]+(sin(xt[ixϕ])-sin(x[ixϕ]))*x[ixp[2]]) + Δt*σβ*n[inβ] #β
  xt[ixr] = x[ixr]                      #r
  xt[ixa] = a                           #a
  return
end

@inbounds function signalUpdateOld!(x, t, n, xt, i)
  a = n[ina]
  xt[ixp] = x[ixp] + Δt*x[ixv] + Δt^2/2*a #p
  xt[ixv] = x[ixv] + Δt*a                 #v
  xt[ixϕ] = x[ixϕ] + Δt*σϕ*n[inϕ]         #ϕ
  xt[ixβ] = x[ixβ] + Δt*σβ*n[inβ]         #β
  xt[ixr] = x[ixr]                        #r
  xt[ixa] = a                             #a
  return
end

@inbounds function signalOut(x, t, n, y)
  signalOutNoNoiseFast(x, t, y)
  y[1:2] += σym*n[iny]
  y[3:4] += σam*n[inya]
  return
end

@inbounds function signalOutRealistic(x, t, n, y)
  signalOutNoNoiseRealistic(x, t, y)
  y[1:2] += σym*n[iny]
  y[3:4] += σam*n[inya]
  return
end

function signalOutNoNoise(x, t, y)
  u(v) = [cos(v) sin(v)]'
  p = [x[1]; x[2]]
  #OBS ADDED MINUS
  complexy = sum(x[ixr].*exp(2π/λ*im*x[ixβ]).*exp(2π/λ*im*(-(p'*u(x[ixϕ]))')))
  y[:] = [real(complexy);
  imag(complexy);
  x[ixa]]
  return
end

@inbounds function signalOutNoNoiseFast(x, t, y)
  realExp = x[ixβ] + cos(x[ixϕ])*x[1] + sin(x[ixϕ])*x[2]
  realExp *= -2π/λ
  y[1] = x[ixr[1]]*cos(realExp[1])+x[ixr[2]]*cos(realExp[2])  #Real part
  y[2] = x[ixr[1]]*sin(realExp[1])+x[ixr[2]]*sin(realExp[2])  #Imag part
  y[3:4] = x[ixa]
  return
end

@inbounds function signalOutNoNoiseRealistic(x, t, y)
  sβ = [norm(s1), norm(s2)]
  #Add phase noise? x[ixβ]
  realExp = - sβ + [norm(s1-x[ixp]), norm(s2-x[ixp])]
  realExp *= 2π/λ
  y[1] = x[ixr[1]]*cos(realExp[1])+x[ixr[2]]*cos(realExp[2])  #Real part
  y[2] = x[ixr[1]]*sin(realExp[1])+x[ixr[2]]*sin(realExp[2])  #Imag part
  y[3:4] = x[ixa]

  #Set the statesto the "right" values
  x[ixϕ] = [atan2(s1[2]-x[ixp[2]], s1[1]-x[ixp[1]]), atan2(s2[2]-x[ixp[2]], s2[1]-x[ixp[1]])]
  x[ixβ] = - sβ + [norm(s1-x[ixp]), norm(s2-x[ixp])]
  return
end

function signalpxy(yhat,y,t)
  w = -(norm([yhat[1],yhat[2]]-[y[1],y[2]]))^2/(2*(σym)^2)#-norm(yhat[3]-y[3])^2/(2*(σam)^2)
end

function sinGen!(u, t)
  randn!(slice(u,3:size(u,1),:))
  A = 10/π
  w = 2π/10
  u[1] = -A/3*w^2*cos(w*Δt*t)
  u[2] = -A*w^2*cos(w*Δt*t/1.5)
  return
end

#Generation of noise for particles given accelerator measurement
function randGen!(u, y, t)
  randn!(u)
  randSet!(u, y, t)
  return
end

function randGenSQMC!(u, y, t)
  randSet!(u, y, t)
  return
end

function randSet!(u, y, t)
  #Method 1
  #u[1,:] *= σam
  #u[1,:] += y[t]
  #Change weights too here

  #Method 2
  u[:,1:2] *= σam*σa/sqrt(σam^2+σa^2)
  u[:,1] += y[3,t]*σa^2/(σam^2+σa^2)
  u[:,2] += y[4,t]*σa^2/(σam^2+σa^2)
  return
end

function twoBeams2DSettings(improvement = true, sinus = true)
  ϕ10 = atan2(s1[2], s1[1])
  ϕ20 = atan2(s2[2], s2[1])
  x0 = [0, 0, 0, 0, ϕ10, ϕ20, 0., 0. , .8, 1.2, 0, 0]
  Ny, NuGen, NuSim = 4, 10, 6
  generator = sinus ? sinGen! : ((u,t) -> rands!(u))
  #est = (x,w,t) -> x[:,findmax(w)[2]]
  est = (x,w,t) -> sum(x.*exp(w'),2)
  f = improvement ? signalUpdate! : signalUpdateOld!
  FilterSettings(f, signalOutNoNoiseFast, signalOutRealistic,
    signalpxy, Ny, generator, (randGen!, :aux), x0, NuGen, NuSim, est)
end

function twoBeams2DSettingsSQMC(;kwargs...)
  s = twoBeams2DSettings(kwargs...)
  s.randsfSim = (randGenSQMC!, :aux, :randn)
  return s
end
