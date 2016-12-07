# some tools
minus = lambda do |f,g|
    lambda do |x|
        f[x] - g[x]
    end
end

div = lambda do |f,g|
    lambda do |x|
        f[x]/g[x]
    end
end

mult = lambda do |f,g|
    lambda do |x|
        f[x] * g[x]
    end
end

norm = lambda do |f|
    lambda do |x|
        f[x].abs
    end
end

const = lambda do |c|
    lambda do |x|
        c
    end
end

plus_eps = lambda do |eps|
    lambda do |f|
        lambda do |x|
            f[x+eps]
        end
    end
end

min_eps = lambda do |eps|
    lambda do |f|
        lambda do |x|
            f[x-eps]
        end
    end
end

inf = lambda do |f,g|
    lambda do |x|
        f[x] < g[x]
    end
end

# The limit operator
lim = lambda do |f,eps,prec|
    lambda do |x|
        if inf[norm[minus[plus_eps[eps/2.0][f], plus_eps[eps][f]]], const[prec]][x]
            plus_eps[eps][f][x]
        else
            lim[f,eps/2.0,prec][x]
        end
    end
end

# The derivative scheme
derivative_sheme = lambda do |f|
    lambda do |x|
        lambda do |eps|
            div[minus[plus_eps[eps][f], min_eps[eps][f]], mult[const[2], const[eps]]][x]
        end
    end
end

# The derivative operator at precision 1e-16 is the limit of the derivative scheme
ddx = lambda do |f|
    lambda do |x|
        lim[derivative_sheme[f][x],1,1e-16][0]
    end
end

f = lambda do |x|
    1000000*Math::E**x + 50000*(Math::E**x - 1)/x - 1726550
end

gNewton = lambda do |x|
    x - f[x]/ddx[f][x]
end

gSecantes = lambda do |x0,x|
    (x0*f[x] - x*f[x0])/(f[x]-f[x0])
end

def getInterval f
    a = 0
    b = 0.2
    until f[a]*f[b] <= 0 || b > 10**7
        a += 0.1
        b += 0.1
    end
    return a, b
end

def getZero g, f
    w0, w = getInterval f

    until (w-w0).abs < 10**-6 || (f[w0]).abs < 10**-6
        w = w0
        w0 = g[w0]
    end

    return w0
end

def getZero2 g, f
    w0, w = getInterval f
    aux = w0

    until (w-w0).abs < 10**-6 || (f[w0]).abs < 10**-6
        w0 = g[w0,w]
        w = aux
        aux = w0
    end

    return w0
end

puts getZero gNewton, f
puts getZero2 gSecantes, f
