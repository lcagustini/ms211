require "matrix"

ddvar = lambda do |f,var|
    lambda do |x,y,z|
        case var
        when "x"
            (f[x+1e-15,y,z] - f[x,y,z])/(1e-15)
        when "y"
            (f[x,y+1e-15,z] - f[x,y,z])/(1e-15)
        else
            (f[x,y,z+1e-15] - f[x,y,z])/(1e-15)
        end
    end
end

f1 = lambda do |x,y,z|
    3*x - Math.cos(y*z) - 0.5
end

f2 = lambda do |x,y,z|
    x**2 - 81*(y+0.1)**2 + Math.sin(z) + 1.06
end

f3 = lambda do |x,y,z|
    Math::E**(-1*x*y) + 20*z + (10*Math::PI - 3)/3
end

f = [f1,f2,f3]

x = [1,1,1]

j = Matrix[[ddvar[f1,"x"][x[0],x[1],x[2]],ddvar[f1,"y"][x[0],x[1],x[2]],ddvar[f1,"z"][x[0],x[1],x[2]]],
           [ddvar[f2,"x"][x[0],x[1],x[2]],ddvar[f2,"y"][x[0],x[1],x[2]],ddvar[f2,"z"][x[0],x[1],x[2]]],
           [ddvar[f3,"x"][x[0],x[1],x[2]],ddvar[f3,"y"][x[0],x[1],x[2]],ddvar[f3,"z"][x[0],x[1],x[2]]]]

l, u = j.lup

def getY l, f
    y = Array.new(3)
    y[0] = - f[0][0,0,0]
    y[1] = - f[1][0,0,0] - l.element(1,0)*y[0]
    y[2] = - f[2][0,0,0] - l.element(2,0)*y[0] - l.element(2,1)*y[1]
    y
end

def getX u, y
    x = Array.new(3)
    x[2] = y[0].fdiv(u.element(2,2))
    x[1] = (y[1] - u.element(1,2)*x[2]).fdiv(u.element(1,1))
    x[0] = (y[2] - u.element(0,2)*x[2] - u.element(0,1)*x[1]).fdiv(u.element(0,0))
    x
end

puts getX(u, getY(l, f))
