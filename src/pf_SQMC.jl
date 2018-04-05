include("hilbertcurve/hilbert_curve.jl")

function resample(s, w)
  # Samples new particles based on their weights. If you find algorithmic optimizations to this rou$
  N = length(w)
  bins = [0.; cumsum(exp(w))]
  j = zeros(Int64,N)
  bo = 1
  for i = 1:N
    for b = bo:N
      #WHY CAN I NOT SKIP THE RIGHT PART!?
      if bins[b] <= s[i] < bins[b+1]
        j[i] = b
        bo = b
        break
      end
    end
  end
  return j
end

function resample(w)
  N = length(w)
  bins = [0.; cumsum(exp(w))]
  s = collect((rand()/N+0):1/N:bins[end])
  resample(s,w)
end

""" `xhat = SQMC(f, g, ginv, x0, y, est, N=200, Nu=lenth(x0))`

Gives estimation of the states `x(:,1),...x(:,T)` given outputs `y(:,1)...y(:,T)` from system:

x(t) = f(x(t-1),t,u) \n
y = gn(x,v),\n
where `g(x,t)=E[gn(x,t,u)]`, `ginv(yhat, y, t) = p(y=gn(x,u)|yhat=g(x))`,
u is uniform noise on `[0,1)`, with `E[x(0)] = xhat0`,
and some estimator `est(xₚ,w,t)` that ouputs the estimate `xhat(:,t)`
given the particles `xₚ` and weights `w`
"""
function SQMC(f!, g!, ginv, randsf, xhat0, y, est; N = 256, Nu = length(xhat0), debug = false, xreal=xreal)
  aux= false
  normal = false
  if isa(randsf, Tuple)
    rands! = randsf[1]
    if randsf[2] == :aux
      aux = true
    end
    if length(randsf)>2
      if randsf[3] != :rand && randsf[3] != :randn
        error("Only allowing :rand or :randn as rands argument")
      elseif randsf[3] == :randn
        normal = true
      end
    end
  else
    rands! = randsf
  end
  if !isinteger(log(2,N))
    error("Number of particles must be a power of 2")
  end
  if debug
    pdata = Void
  end
  T = size(y,2)
  #Estimates at time 1:T
  xhat = Array{Float64,2}(length(xhat0), T)
  #Current predicted output for each particle
  yhat = Array{Float64,2}(size(y,1), N)
  x = Array{Float64,2}(length(xhat0), N)
  hilbertIndex = Array{Int64}(size(x,2))
  xtemp = similar(x)
  w = Array{Float64,1}(N)
  #x[:,1] = xhat0; #Hack to use the pfStep fuction. Works by letting `a` be ones
  x = repmat(xhat0,1,N);

  #Step (a) (Draw u and x), OBS u will be of size NxNu
  u, = sobol(Nu,N)
  if normal
    unit2Normal!(u)
  end
  if !aux
    rands!(u)
  else
    rands!(u,y,1)
  end

  if debug
    p = Void
    xOld = similar(x)
    xOld[:] = x[:]
  end
  pfStep!(f!,g!,ginv,x,xtemp,u,ones(Int,N),1:N,1,w,yhat,y,N, 1:Nu)
  xhat[:,1] = est(x,w,1)
  if debug
    #p = plotPoints(x, w, y, yhat, N, ones(Int,N), 1, xreal, xhat, xOld, p)
    p = pploti(x, exp(w)./sum(exp(w)), y, yhat, N, ones(Int,N), 1, xreal, xhat, xOld, p, density = true,leftOnly=true, xIndices = [1], yIndices = [],  slidef=1.0)
    xOld[:] = x[:]
  end
  #u = Array{Float64,2}(N,Nu+1)
  #nextseed, MeM = sobol!(u,Nu+1,N)
  u, nextseed, MeM = sobol(Nu+1,N)
  #The first column in u is for resample purpose
  Nus = 2:(Nu+1)
  for t = 2:T
    #nextseed = sobol!(u, nextseed, MeM)
    nextseed = sobol!(u, nextseed, MeM)
    if normal
      unit2Normal!(slice(u,:,Nus))
    end
    if !aux
      rands!(slice(u,:,Nus))
    else
      rands!(slice(u,:,Nus),y,t)
    end
    τ = sortperm(slice(u,:,1))
    σ = sortperm(slice(ψ!(x,hilbertIndex),:))
    a = resample(slice(u,τ,1), slice(w,σ))
    # Time update, no sigma here?
    #if randsf[4]
    #Make this faster
      x = x[:,σ]
      w = w[σ]
    #end
    pfStep!(f!,g!,ginv,x,xtemp,u,a,τ,t,w,yhat,y,N, Nus)
    xhat[:,t] = est(x,w,1)
    if debug
      #p = plotPoints(x, w, y, yhat, N, a, t, xreal, xhat, xOld, p)
      p = pploti(x, exp(w)./sum(exp(w)), y, yhat, N, a, t, xreal, xhat, xOld, p, density = true,leftOnly=true, xIndices = [1], yIndices = [],  slidef=1.0)
      xOld[:] = x[:]
    end
  end
  xhat
end

function unit2Normal!(u)
  for j = 1:size(u,2)
    for i = 1:size(u,1)
      u[i,j] = sqrt(2)*erfinv(2*u[i,j]-1);
    end
  end
end

#Should only be used for when 0<x<1
function ψUnit!(x,hbtIdx)
  #Figure out a resonable resolution (Number of indices has to be lower than max int)
  ncoords = size(x,1)
  N=floor(Int64,2^(60/ncoords))-1 #The resolution in each direction
  #Vector to keep the integer (scaled) value for each coordinate
  xInts = Array{Int64}(size(x,1))
  for j = 1:size(x,2)
    for i = 1:size(x,1)
      xInts[i] = fld(x[i,j]*N,1/N)
    end
    hbtIdx[j] = hilbert2int(xInts)
  end
  hbtIdx
end

function ψ!(x,hbtIdx)
  #Figure out a resonable resolution (Number of indices has to be lower than max int)
  ncoords = size(x,1)
  N=floor(Int64,2^(60/ncoords))-1 #The resolution in each direction
  #Find a reasonable scale
  maxVal = maximum(x,2)
  minVal = minimum(x,2)
  #Only use incices for which max!=min in sorting
  idx = find((v)->v!=0,maxVal-minVal)
  length(idx) == 0 && error("Error in sortnig, all values equal")
  #Vector to keep the integer (scaled) value for each coordinate
  xInts = Array{Int64}(length(idx))
  for j = 1:size(x,2)
    for (i,val) in enumerate(idx)
      xInts[i] = fld((x[val,j]-minVal[val])*N,maxVal[val]-minVal[val])
    end
    hbtIdx[j] = hilbert2int(xInts)
  end
  hbtIdx
end

#Randsf is either a function that fills the u matrix (only argument) with random numebrs or a tuple with a function and a symbol.
#If symbol is :aux then the function should accept the matrix as well as the measurement and time as arguments
function pf(f!, g!, ginv, randsf, xhat0, y, est; N = 200, Nu = length(xhat0), debug=false, xreal=xreal)
  aux= false
  if isa(randsf, Tuple)
    rands! = randsf[1]
    if randsf[2] == :aux
      aux = true
    end
  else
    rands! = randsf
  end
  T = size(y,2)

  #Estimates at time 1:T
  xhat = Array{Float64,2}(length(xhat0), T)
  #Current predicted output for each particle
  yhat = Array{Float64,2}(size(y,1), N)
  x = Array{Float64,2}(length(xhat0), N)
  xtemp = similar(x)
  w = Array{Float64,1}(N)
  τ = 1:N
  x = repmat(xhat0,1,N); #Hack to use the pfStep fuction. Works by letting `a` be ones

  u = Array{Float64,2}(N,Nu)
  if !aux
    rands!(u)
  else
    rands!(u,y,1)
  end

  if debug
    p = Void
    xOld = similar(x)
    xOld[:] = x[:]
  end
  pfStep!(f!,g!,ginv,x,xtemp,u,ones(Int,N),τ,1,w,yhat,y,N, 1:Nu)
  xhat[:,1] = est(x,w,1)
  if debug
    #p = plotPoints(x, w, y, yhat, N, ones(Int,N), 1, xreal, xhat, xOld, p)
    p = pploti(x, exp(w)./sum(exp(w)), y, yhat, N, ones(Int,N), 1, xreal, xhat, xOld, p, density = true, leftOnly=true, xIndices = [1], yIndices = [], slidef=1.0)
    xOld[:] = x[:]
  end
  for t = 2:T
    #Avoiding SQMC resample because of gaussian noise (change back?)
    a = resample(w)
    #a = resample(slice(u,τ,1), slice(w,:))
    #Rand
    if !aux
      rands!(u)
    else
      rands!(u,y,t)
    end
    # Time update, no sigma here?
    pfStep!(f!,g!,ginv,x,xtemp,u,a,τ,t,w,yhat,y,N,1:Nu)
    xhat[:,t] = est(x,w,1)
    if debug
      #p = plotPoints(x, w, y, yhat, N, a, t, xreal, xhat, xOld, p)
      p = pploti(x, exp(w)./sum(exp(w)), y, yhat, N, a, t, xreal, xhat, xOld, p, density = true,leftOnly=true, xIndices = [1], yIndices = [],  slidef=1.0)
      xOld[:] = x[:]
    end
  end
  xhat
end

#Time update and weighting. Chances `x`, `yhat` and `w`
function pfStep!(f!,g!,ginv,x,xtemp,u,a,τ,t,w,yhat,y,N,uRange)
  for i = 1:N
    xS = slice(x, :, a[i])
    uS = slice(u, τ[i], uRange)
    xtS = slice(xtemp,:,i)
    f!(xS, t, uS, xtS, i)
    g!(xtS, t, slice(yhat,:,i))
    #println(ginv(slice(yhat,:,i), slice(y,:,t), t))
    w[i] = ginv(slice(yhat,:,i), slice(y,:,t), t)
  end
  x[:] = xtemp[:]
  offset = maximum(w)
  mySum = sum(exp(w-offset))
  normConstant = log(mySum)+offset
  w[:] = w[:] - normConstant
end

function generateRealSequence(f!, gn!, x0, T, rands!, Ny, Nu)
  x = Array{Float64,2}(length(x0), T)
  xt = Array{Float64,2}(length(x0), 1)
  u = Array{Float64,1}(Nu)
  y = Array{Float64,2}(Ny, T)
  rands!(u, 1)
  f!(x0, 1, u,xt,1)
  x[:,1] = xt[:]
  gn!(slice(x,:,1), 1, u, slice(y,:,1))
  for t = 2:T
    rands!(u, t)
    f!(x[:,t-1], t, u,xt,t)
    x[:,t] = xt[:]
    gn!(slice(x,:,t), t, u, slice(y,:,t))
  end
  x, y
end

type FilterSettings
  f::Function
  g::Function
  gn::Function
  ginv::Function
  Ny::Int
  randsfGen
  randsfSim
  x0::Vector
  NuGen::Int
  NuSim::Int
  est
end

function runTest(N, method, s, T, debug = false)
  x, y = generateRealSequence(s.f, s.gn, s.x0, T, s.randsfGen, s.Ny, s.NuGen)
  xhat = method(s.f, s.g, s.ginv, s.randsfSim, s.x0, y, s.est,
  N = N, Nu = s.NuSim, debug = debug, xreal = x)
  x, xhat
end
