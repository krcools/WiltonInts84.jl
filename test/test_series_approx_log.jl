@testitem "q/a small q/b large" begin
    a = -0.01927176760943159
    b = 0-1.834267662268588e-6
    p = 2.7959608455489207e-6
    h = 0.0
    m = [0.008954817714265238, 0.0, -0.9999599048160402]
    UB = Val{1}

    T = typeof(a)
    a2 = a^2
    b2 = b^2
    q2 = p^2 + h^2

    ra = sqrt(a2 + q2)
    rb = sqrt(b2 + q2)

    WiltonInts84.segintsg(a, b, p, h, m, UB)
    j1 = log(a/b) + log((1-(q2/b2)/4) / (1-(q2/a2)/4))
    j2 = log(b + rb) - log(a + ra)
    j3 = log((b + rb) / (a * (-0.5*q2/a2 + 0.125*(q2/a2)^2)))

    @test abs(j3-j1) > 1e-1
    @test abs(j3-j2) < sqrt(eps(T))
end