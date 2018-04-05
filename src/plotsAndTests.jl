include("../particlePlot/particlePlot.jl")

# State x = [pₓ, vₓ , ϕ₁, ϕ₂, β₁, β₂, r₁, r₂, a]
# Noise u = [a, ϕ₁, ϕ₂, β₁, β₂]
# Out   y₁ = ∑rₖe^(iβₖ)e^(i2π<pₖ,u(ϕₖ)>)+ny
#       y₂ = a + na
#
# function runTest(N, method, debug, T = 20, withImprovement=true)
#   if withImprovement
#     f! = signalUpdate!
#   else
#     f! = signalUpdateOld!
#   end
#   #g = signalOutNoNoise
#   g = signalOutNoNoiseFast
#   gn = signalOut
#   Ny = 3
#   ginv = signalpxy
#   randsf = (randGen!, :aux)
#   #est = (x,w,t) -> x[:,findmax(w)[2]]
#   est = (x,w,t) -> sum(x.*exp(w'),2)
#
#   x0 = [0., 0, .2, .7, 0., 0. , .8, 1.2, 0]
#   x, y = generateRealSequence(f!, gn, x0, T, randn!, Ny, Nu = 8)
#   xhat = method(f!, g, ginv, randsf, x0, y, est, N = N, Nu = 5, debug=debug, xreal=x)
#   x, xhat
# end
#
function rms(x)
  sqrt(1/length(x)*sum(x.^2))
end

function testRadio()
  debug = true
  Ns = 2.^(5:5)
  M = 1
  T = 100
  RMS = Array{Float64}(length(Ns),M,2)
  largeRMS = Array{Float64}(length(Ns),2)
  rmsMean = Array{Float64,2}(length(Ns),2)
  rmsVariance = Array{Float64,2}(length(Ns),2)
  # Method, Particles, Run nbr, variable, time
  allxhats = Array{Float64}(2,length(Ns),M,4,T)
  #for (methodidx,method) in enumerate([false,true])
  for (methodidx,method) in enumerate([pf,SQMC])
    println("Method: "*string(method))
    xhat, xreal = 0, 0
    @time for (i, N) in enumerate(Ns)
      println(string(N)*", ")
      rmslocal = Array{Float64,1}(M)
      for j = 1:M
        print(string(j)*",")
        #xreal, xhat = runTest(N,pf,debug,T,method)
        #xreal, xhat = runTest(N,pf,twoBeams2DSettings(method),T,debug)
        if method == pf
          xreal, xhat = runTest(N,pf,twoBeams2DSettings(),T,debug)
        else
          xreal, xhat = runTest(N,SQMC,twoBeams2DSettingsSQMC(),T,debug)
        end
        allxhats[methodidx,i,j,1:4,:] = xhat[[1,2,5,6],:]
        RMS[i,j,methodidx] = rms(xreal[1,:]-xhat[1,:])
      end
    end
    rmsMean[:,methodidx] = mean(RMS[:,:,methodidx],2)
    rmsVariance[:,methodidx] = std(RMS[:,:,methodidx],2)./sqrt(M)
    for (i, N) in enumerate(Ns)
      largeRMS[i,methodidx] = length(find(RMS[i,:,methodidx].>rmsMean[i,methodidx]+2*rmsVariance[i,methodidx]))
    end
  end
  rmsMean, rmsVariance, Ns, largeRMS, RMS, allxhats
end

function testPlot(rmsMean, rmsVariance, Ns)
  p = plot(x=Ns,y=rmsMean[:,1],c=:red)
  plot!(x=Ns,y=rmsMean[:,2],c=:blue)
  plot!(x=Ns,y=rmsMean[:,1]+rmsVariance[:,1]*2,c=:red)
  plot!(x=Ns,y=rmsMean[:,1]-rmsVariance[:,1]*2,c=:red)
  plot!(x=Ns,y=rmsMean[:,2]+rmsVariance[:,2]*2,c=:blue)
  plot!(x=Ns,y=rmsMean[:,2]-rmsVariance[:,2]*2,c=:blue)
  gui(p)
end
# 
# function heatBoxPlot(plt, x, t, minmax, nbinsy=30)
#   if maximum(x)-minimum(x) > 0
#     heatmap!(plt, x=(t-1.5):1/length(x):(t-0.5),y=x',nbins=(1,nbinsy), ylims=tuple(minmax...))
#   end
#   return
# end
# 
# function distPlotY(plt, x, t, minmax, color=:blue, nbinsy=30)
#   if maximum(x)-minimum(x) > 0
#     y, nbr = hist(x', nbinsy)
#     xval = nbr[:]#cumsum(nbr[:])
#     plot!(plt, [0, xval[:]]/xval[end], y, linecolor=color, ylims=tuple(minmax...))
#   end
#   return
# end
# 
# """ `plotPoints(x, w, y, yhat, N, a, t, xreal, xhat, xOld, pdata)`
# 
# To be called inside a particle filter, plots either particle density (`density=true`) or individual particles (`density=false`) \n
# To only plot the particle trajectories, set (`leftOnly=false`)\n
# Will plot all the real states in `pltxIdx` as well as the expected vs real measurements of `pltyIdx`.
# 
# Arguments: \n
#   * `x`: `Array(M,N)`. The states for each patyicle where `M` number of states, `N` number of Particles
#   * `w`: `Array(N)`. logprobability of each particle
#   * `y`: `Array(R,T)`. All true outputs. `R` is number of outputs, `T` is total number of time steps (will only use index `t`)
#   * `yhat`: `Array(R,N)` The expected output per particle. `R` is number of outputs, `N` number of Particles
#   * `N`, Number of particles
#   * `a`, `Array(N)`, reorderng of particles (e.g. `1:N`)
#   * `t`, Current time step
#   * `xreal`: `Array(M,T)`. All true states. `R` is number of states, `T` is total number of time steps (will only use index `t`)
#   * `xhat`: Not used
#   * `xOld`: Same as `x`, but for previous time step, only used when `!density` to show states origins
#   * `pdata`: Persistant data for plotting. Set to void in first call and pdataOut on remaining \n
# Returns: `pdataOut`
# """
# function plotPoints(x, w, y, yhat, N, a, t, xreal, xhat, xOld, pdata)
#   immerse()
#   leftOnly = false
#   density = true
#   cols = leftOnly?1:2
#   grd = (r,c) -> (r-1)*cols+c
#   println("Surviving: "*string((N-length(setdiff(Set(1:N),Set(a))))/N))
#   vals = [x;yhat]
#   realVals = [xreal;y]
#   valsOld = [xOld;yhat]
#   pltxIdx = [1,2,5, 6]# 7 9]
#   pltyIdx = [1,2]
#   pltIdx = [pltxIdx; size(x,1)+pltyIdx]
#   println(pltIdx)
#   if pdata == Void
#     pdata = (subplot(layout=cols*ones(Int,length(pltIdx))), Array{Float64,2}(length(pltIdx),2))
#     gui(pdata[1])
#   end
#   p, minmax = pdata
#   println(size(vals[pltIdx,:],2))
#   minmax = [min(minmax[:,1], minimum(vals[pltIdx,:],2)) max(minmax[:,2], maximum(vals[pltIdx,:],2))]
#   #c = exp(w[:]-minimum(w))*3
#   c = exp(w[:])*5*N
#   for (i, val) in enumerate(pltIdx)
#     if !leftOnly
#       oldFig = p.plts[grd(i,2)].o[1]
#       newPlot = plot()
#     end
#     if density
#       #Plot the heatmap on the left plot
#       heatBoxPlot(p.plts[grd(i,1)], vals[val,:], t, minmax[i,:])
#     end
#     for j = 1:N
#       if !density
#         #Plot the line on the left plot
#         plot!(p.plts[grd(i,1)], [t-1.5,t-1], [valsOld[val,j], valsOld[val,j]], legend=false)
#         plot!(p.plts[grd(i,1)], [t-1,t-0.5], [valsOld[val,a[j]], vals[val,j]], legend=false)
#       end
#       if !leftOnly && !density
#         #Plot each of the dots on the right side
#         plot!(newPlot, [j,j], [vals[val,j]',vals[val,j]'], marker=(:circle, :red, c[j])  , legend=false)
#       end
#     end
#     if !leftOnly
#       if density
#         distPlotY(newPlot, vals[val,:] , t, minmax[i,:], :blue)
#         #Weight can not be multiplied with value! DUH!
#         #distPlotY(newPlot, vals[val,:].*w' , t, :red)
#       end
#       #Fix the bug with updating plot
#       p.plts[grd(i,2)] = copy(newPlot)
#       p.plts[grd(i,2)].o = (oldFig, p.plts[grd(i,2)].o[2])
#     end
#     #Plot Real State Here
#     plot!(p.plts[grd(i,1)], (t-1):t, [realVals[val,t], realVals[val,t]], legend=false, color=:black, linewidth=5)
#   end
#   gui(p)
#   print("Waiting for command. q to Quit, ^D to run all:\n")
#   if readline(STDIN) == "q\n"
#     error("Exit command")
#   end
#   (p, minmax)
# end
