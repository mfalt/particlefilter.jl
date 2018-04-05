using Gadfly
using Colors

include("pf_SQMC.jl")
include("sobol.jl")


function GordonKitagawaUpdate!(x,t, u, xtemp, i)
  b₁, b₂, b₃, b₄, σ = .5, 25, 8, 1.2, sqrt(1);
  xtemp[1,i] = b₁*x[1]+b₂*x[1]/(t+x[1]^2)+b₃*cos(b₄*t)+σ*erfinv(u[1]*2-1);
  return
end

function GordonKitagawaOut(x,t)
  σ = sqrt(.1);
  a = 20
  x[1]^2/a+σ*randn()[1]
end

function GordonKitagawaOutNoNoise(x, t)
  a = 20
  x[1].^2./a
end

function GordonKitagawapxy(yhat,y,t)
  σ = sqrt(.1);
  #w = 1/sqrt(2π)*e^(-(yhat[1]-y[1])^2/2)
  w = -(yhat[1]-y[1])^2/(2*σ^2)
end



function generateRealSequence(f!, g, x0, T, N = length(x0))
  x = Array{Float64,2}(length(x0), T)
  f!(x0, 1, rand(N),x,1)
  y1 = g(x[:,1],1)
  y = Array{Float64,2}(length(y1), T)
  y[:,1] = y1
  for t = 2:T
    f!(x[t-1], t, rand(N),x,t)
    y[:,t] = g(x[t], t)
  end
  x, y
end

function estMean(x,w,t)
  sum(x[:].*exp(w[:]))
end

function runTest(N, method, debug)
  T = 150

  ##Gordon Kitagawa
  f! = GordonKitagawaUpdate!#(x,t,u,xtemp,i) -> GordonKitagawaUpdate!(x[1],t,u,xtemp,i)
  g = GordonKitagawaOutNoNoise#(x,t) -> GordonKitagawaOutNoNoise(x[1])
  gn = GordonKitagawaOut#(x,t) -> GordonKitagawaOut(x[1])
  ginv = GordonKitagawapxy#(yhat, y,t) -> GordonKitagawapxy(yhat,y)
  rands = rand
  ##LTI System
  #σ = .5
  #f! = (x,t,u,xtemp,i) -> xtemp[1,i] = .8x[1]+4*erfinv(u[1]*2-1)
  #g = (x,t) -> 2*x[1]
  #gn = (x,t) -> 2*x[1] + σ*randn()[1]
  #ginv = (yhat, y, t) -> -(yhat[1]-y[1])^2/(2*σ^2)

  ##Estimator
  #est = (x,w,t) -> x[findmax(w)[2]]
  est = (x,w,t) -> sum(x[:].*exp(w[:]))

  x0 = .5
  x, y = generateRealSequence(f!, gn, x0, T)
  xhat = method(f!, g, ginv, rands, x0, y, est, N, debug=debug, xreal=x)
  x, xhat
end



function plotPoints(x, w, y, N, a, τ, t, xreal, xhat)
  c = w[:]-minimum(w)+1
  ##Use for GordonK
  #p = plot(layer(x=collect(1:N), y=x[:], Geom.point, color=c),
  #        layer(x=[1,N],y=ones(2).*sqrt(20*max(y[:,t],0)),Geom.line, Theme(default_color=color(colorant"red"))),
  #        layer(x=[1,N],y=-ones(2).*sqrt(20*max(y[:,t],0)),Geom.line, Theme(default_color=color(colorant"red"))),
  #        layer(x=[1,N],y=ones(2).*xreal[:,t],Geom.line, Theme(default_color=color(colorant"blue"),line_width=2px)),
  #        layer(x=[1,N],y=ones(2).*xhat[:,t],Geom.line, Theme(default_color=color(colorant"black"),line_width=4px)),
  #        Guide.XLabel("Particle "*string(t)), Guide.YLabel("Estimate"), Coord.Cartesian(ymin=-15,ymax=15))
  ##Use for LTI
  p = plot(layer(x=collect(1:N), y=x[:], Geom.point, color=c),
  layer(x=[1,N],y=ones(2)*1/2.*y[:,t],Geom.line, Theme(default_color=color(colorant"red"))),
  layer(x=[1,N],y=ones(2).*xreal[:,t],Geom.line, Theme(default_color=color(colorant"blue"),line_width=2px)),
  layer(x=[1,N],y=ones(2).*xhat[:,t],Geom.line, Theme(default_color=color(colorant"black"),line_width=4px)),
  Guide.XLabel("Particle "*string(t)), Guide.YLabel("Estimate"), Coord.Cartesian(ymin=-10,ymax=10))

  display(p)
  print("here")
  readline(STDIN)
end


function rms(x)
  sqrt(1/length(x)*sum(x.^2))
end

function testSQMC()
  debug = false
  Ns = 2.^(2:7)
  M = 50
  RMS = Array{Float64,2}(length(Ns),M)
  largeRMS = Array{Float64}(length(Ns),2)
  rmsMean = Array{Float64,2}(length(Ns),2)
  rmsVariance = Array{Float64,2}(length(Ns),2)
  for (methodidx,method) in enumerate([SQMC,pf])
    xhat, xreal = 0, 0
    @time for (i, N) in enumerate(Ns)
      rmslocal = Array{Float64,1}(M)
      for j = 1:M
        xreal, xhat = runTest(N,method,debug)
        RMS[i,j] = rms(xreal-xhat)
      end
    end
    rmsMean[:,methodidx] = mean(RMS,2)
    rmsVariance[:,methodidx] = std(RMS,2)./sqrt(M)
    for (i, N) in enumerate(Ns)
      largeRMS[i,methodidx] = length(find(RMS[i,:].>rmsMean[i,methodidx]+2*rmsVariance[i,methodidx]))
    end
  end
  rmsMean, rmsVariance, Ns, largeRMS, RMS
end

function testPlot(rmsMean, rmsVariance, Ns)
  p = plot(
  layer(x=Ns,y=rmsMean[:,1],Geom.line,Theme(default_color=color(colorant"red"))),
  layer(x=Ns,y=rmsMean[:,2],Geom.line,Theme(default_color=color(colorant"blue"))),
  layer(x=Ns,y=rmsMean[:,1]+rmsVariance[:,1]*2,Geom.line,Theme(default_color=color(colorant"red"))),
  layer(x=Ns,y=rmsMean[:,1]-rmsVariance[:,1]*2,Geom.line,Theme(default_color=color(colorant"red"))),
  layer(x=Ns,y=rmsMean[:,2]+rmsVariance[:,2]*2,Geom.line,Theme(default_color=color(colorant"blue"))),
  layer(x=Ns,y=rmsMean[:,2]-rmsVariance[:,2]*2,Geom.line,Theme(default_color=color(colorant"blue"))),
  )
end


function testPlot(N)
  for j = 1:N
    M = 1000
    vals = Array{Float64,2}(M,2)
    ds = linspace(rand()/M,1,M)[1:end]
    for i = 1:M
      o =  d2xy(ds[i])
      vals[i,:] = [o[1], o[2]]
    end
    plot(vals[:,1], vals[:,2])
  end
end

function testPlot2(N)
  for j = 1:N
    M = 1000
    vals = Array{Float64,2}(M,2)
    ds = sort(rand(M))
    for i = 1:M
      o =  d2xy(ds[i])
      vals[i,:] = [o[1], o[2]]
    end
    plot(vals[:,1], vals[:,2])
  end
end
