include("pf_SQMC.jl")
include("plotsAndTests.jl")
include("sobol.jl")

const Δt = 0.02
#The size of the random acceleration
const σa = 10.0 #m/s^2
#The size of the random change in ϕ
const σϕ = 0.5
#The size of the random change in β
const σβ = 0.05
#The size of the noise on acceleration measurements
const σam = .1
const σym = .1
#Wavelength in meters
const λ = 0.13

# State x = [pₓ, vₓ , ϕ₁, ϕ₂, β₁, β₂, r₁, r₂, a]
# Noise u = [a, ϕ₁, ϕ₂, β₁, β₂, ny1, ny2, na]
# Out   y₁ = ∑rₖe^(iβₖ)e^(i2π<pₖ,u(ϕₖ)>)+ny
#       y₂ = a + na

#OBS LET NOISE BE NORMAL
@inbounds function signalUpdate!(x, t, n, xt, i)
  a = n[1]
  py = 0
  u(v) = [-sin(v) cos(v)]'
  xt[1] = x[1] + Δt*x[2] + Δt^2/2*a     #p⁺ = p + Δt*v + Δt²/2*a
  xt[2] = x[2] + Δt*a                   #v⁺ = v + Δt*a
  #xt[3] = x[3] + Δt*n[2]                  #ϕ₁⁺= ϕ₁ + Δt*nϕ₁
  #xt[4] = x[4] + Δt*n[3]                  #ϕ₂⁺= ϕ₂ + Δt*nϕ₂
  xt[3:4] = x[3:4] + Δt*σϕ*n[2:3]
  #xt[5] = x[5] - (u(xt[3])-u(x[3]))*[x[1]; py] #β₁⁺= β₁ - <p,u(ϕ₁⁺)-u(ϕ₁)> + nβ₁
  #xt[6] = x[6] - (u(xt[4])-u(x[4]))*[x[1]; py] #β₂⁺= β₂ - <p,u(ϕ₂⁺)-u(ϕ₂)> + nβ₂
  #xt[5:6] = x[5:6] - (u(xt[3:4])-u(x[3:4]))*[x[1]; py] + Δt*σβ*n[4:5]
  xt[5:6] = x[5:6] - (-sin(xt[3:4])+sin(x[3:4]))*x[1] + Δt*σβ*n[4:5]
  xt[7:8] = x[7:8]
  xt[9] = a
  return
end

@inbounds function signalUpdateOld!(x, t, n, xt, i)
  a = n[1]
  py = 0
  xt[1] = x[1] + Δt*x[2] + Δt^2/2*a
  xt[2] = x[2] + Δt*n[1]
  xt[3:4] = x[3:4] + Δt*σϕ*n[2:3]
  xt[5:6] = x[5:6] + Δt*σβ*n[4:5]
  xt[7:8] = x[7:8]
  xt[9] = a
  return
end

@inbounds function signalOut(x, t, n, y)
  signalOutNoNoise(x, t, y)
  y[1] += σym*n[6]
  y[2] += σym*n[7]
  y[3] += σam*n[8]
  return
end

function signalOutNoNoise(x, t, y)
  ϕi, βi, ri, ai = 3:4, 5:6, 7:8, 9
  u(v) = [-sin(v) cos(v)]'
  p = [x[1]; 0]
  complexy = sum(x[ri].*exp(2π/λ*im*x[βi]).*exp(2π/λ*im*(p'*u(x[ϕi]))'))
  y[:] = [real(complexy);
  imag(complexy);
  x[ai]]
  return
end

@inbounds function signalOutNoNoiseFast(x, t, y)
  realExp = x[5:6] - sin(x[3:4])*x[1]
  realExp *= 2π/λ
  y[1] = x[7]*cos(realExp[1])+x[8]*cos(realExp[2])
  y[2] = x[7]*sin(realExp[1])+x[8]*sin(realExp[2])
  y[3] = x[9]
  return
end

function signalpxy(yhat,y,t)
  #w = 1/sqrt(2π)*e^(-(yhat[1]-y[1])^2/2)
  #Without Aux
  #w = -(norm([yhat[1];yhat[2]]-[y[1];y[2]]))^2/(2*(σym)^2)-norm(yhat[3]-y[3])^2/(2*(σam)^2)
  #With Aux
  w = -(norm([yhat[1],yhat[2]]-[y[1],y[2]]))^2/(2*(σym)^2)#-norm(yhat[3]-y[3])^2/(2*(σam)^2)
end

function sinGen!(u, t)
  randn!(slice(u,2:size(u,1),:))
  A = 10/π
  w = 2π/10
  u[1] = -A*w^2*cos(w*Δt*t)
  return
end

#Generation of noise for particles given accelerator measurement
function randGen!(u, y, t)
  randn!(u)
  #Method 1
  #u[1,:] *= σam
  #u[1,:] += y[t]
  #Change weights too here

  #Method 2
  u[:,1] *= σam*σa/sqrt(σam^2+σa^2)
  u[:,1] += y[3,t]*σa^2/(σam^2+σa^2)
  return
end

function twoBeamsOnLineSettings(improvement = true, sinus = true)
  x0 = [0., 0, .2, .7, 0., 0. , .8, 1.2, 0]
  Ny, NuGen, NuSim = 3, 8, 5
  generator = sinus ? sinGen! : ((u,t) -> rands!(u))
  #est = (x,w,t) -> x[:,findmax(w)[2]]
  est = (x,w,t) -> sum(x.*exp(w'),2)
  f = improvement ? signalUpdate! : signalUpdateOld!
  FilterSettings(f, signalOutNoNoiseFast, signalOut,
    signalpxy, Ny, generator, (randGen!, :aux), x0, NuGen, NuSim, est)
end
