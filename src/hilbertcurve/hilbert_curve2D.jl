#Convert number in [0,1)^2 to [0,1)
function xy2d(x::Float64, y:: Float64)
    N = 2^30
    xint = floor(Int64, x*N)
    yint = floor(Int64, y*N)
    # d will be in [0,N^2)
    d = xy2d(N, xint, yint)/N^2
end

function d2xy(d:: Float64)
    N::Int64 = 2^30
    dint = floor(Int64, d*N^2)
    x::Float64, y::Float64 = d2xy(N^2, dint)
    x/N, y/N
end

#Converted from wikipedia
#convert (x,y) to d
function xy2d(n::Integer, x::Integer, y::Integer)
    rx, ry, s, d = 0, 0, 0, 0
    s = div(n,2)
    while s>0
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * ((3*rx) $ ry)
        x, y = rot(s, x, y, rx, ry)
        s = div(s,2)
    end
    return d::Integer
end

#convert d to (x,y)
function d2xy(n::Integer, d::Integer)
    rx, ry, s, t = 0, 0, 0, d
    x, y = 0, 0
    s = 1
    while s<n
        rx = 1 & (div(t,2))
        ry = 1 & (t $ rx)
        x, y = rot(s, x, y, rx, ry)
        x += s*rx
        y += s*ry
        t = div(t,4)
        s *= 2
    end
    return x::Integer, y::Integer
end

#rotate/flip a quadrant appropriately
function rot(n::Integer, x::Integer, y::Integer, rx::Integer, ry::Integer)
    if ry == 0
        if rx == 1
            x = n-1 - x
            y = n-1 - y
        end
        #Swap x and y
        t = x
        x = y
        y = t
    end
    return x::Integer, y::Integer
end
