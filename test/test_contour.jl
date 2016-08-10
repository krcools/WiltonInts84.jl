using Base.Test
using WiltonInts84
using FixedSizeArrays

WI = WiltonInts84

p1 = Vec(1.0, 0.0, 0.0)
p2 = Vec(0.0, 1.0, 0.0)
p3 = Vec(0.0, 0.0, 0.0)

c = (p1+p2+p3)/3
r, R = 0.35, 0.45

q = WI.contour(p1,p2,p3,c,r,R)

@test length(q.segments) == 6
@test length(q.arcs) == 6
@test length(q.circles) == 0
