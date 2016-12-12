# MS211 - Exerc 3
# Lucas Cunha Agustini - 172655
# Andre Figueiredo de Almeida - 164047
# Joao Víctor Brazileu Spuri - 155943

a = [[32, 8, 8, 0, 0, 0, 0, 0, 0, 0],
     [4, 25, -3, 6, 0, 0, 0, 10, 0, 0],
     [4, -3, 30, 8, 0, 2, 0, -9, 0, 0],
     [0, 3, 4, 23.5, 4, -1, 0, 0, 0, 0],
     [0, 0, 0, 8, 17, 3, 2, -1, 0, 0],
     [0, 0, 2, -2, 3, 15, 1, -3, 3, 0],
     [0, 0, 0, 0, 2, 1, 21, -4, -7, 0],
     [0, 10, -9, 0, -1, -3, -4, 52, -1, 4],
     [0, 0, 0, 0, 0, 3, -7, -1, 14, 0],
     [0, 0, 0, 0, 0, 0, 0, 4, 0, 24]]

b = [28, -24.5, 87.7, 71.25, 7.9, 47.35, -46.45, 8.45, -3.5, 80]

c = Array.new(10) { Array.new(10) }

g = Array.new(10)

for i in (0...10) do
    for j in (0...10) do
        if i == j
            c[i][j] = 0
        else
            c[i][j] = -a[i][j].fdiv(a[i][i])
        end
    end
end

for i in (0...10) do
    g[i] = b[i].fdiv(a[i][i])
end

x = Array.new(10) {10}
y = Array.new(10) {10}

def precise? x, x0, a, b
    p = Array.new(x.length)
    for i in (0...x.length) do
        p[i] = (x[i] - x0[i]).abs
    end
    p.max.fdiv(x.max.abs) < 10**-8
end

# Gauss-Jacobi
# Utiliza a matriz original para calcular k
(10**2).times do
    l = x.clone
    for i in (0...x.length) do
        k = 0
        for j in (0...x.length) do
            k += c[i][j]*l[j]
        end
        x[i] = g[i] - k
    end
    break if precise? x, l, a, b
end

# Gauss-Seidel
# Utiliza a matriz que ja esta sendo gerada para calcular k
(10**2).times do
    l = x.clone
    for i in (0...x.length) do
        k = 0
        for j in (0...x.length) do
            k += c[i][j]*y[j]
        end
        y[i] = g[i] - k
    end
    break if precise? y, l, a, b
end

print x
puts ""
print y
