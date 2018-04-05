include("twoBeams2D.jl")
N = 20
xhats = Array{Float64}(12,1000,N,2)
x = 0
first = true
p1 = plot(title="Improved, N = 3000")
p2 = plot(title="Improved, N = 3000")
for i = 1:N
  x, xhat = runTest(3000, pf, twoBeams2DSettings(true), 1000)
  xhats[:,:,i,1] = xhat
  if first
    plot!(p1, x[1,:]', x[2,:]', width=3)
    plot!(p2, x[5:6,:]', width=3)
  end
  first = false
  plot!(p1, xhat[1,:]', xhat[2,:]')
  plot!(p2, xhat[5:6,:]')
  gui(p1)
  gui(p2)
  println(i)
end

first = true
p3 = plot(title="Old, N = 3000")
p4 = plot(title="Old, N = 3000")
for i = 1:N
  x, xhat = runTest(4096, pf, twoBeams2DSettings(false), 1000)
  xhats[:,:,i,2] = xhat
  if first
    plot!(p3, x[1,:]', x[2,:]', width=3)
    plot!(p4, x[5:6,:]', width=3)
  end
  first = false
  plot!(p3, xhat[1,:]', xhat[2,:]')
  plot!(p4, xhat[5:6,:]')
  gui(p3)
  gui(p4)
  println(i)
end
