using Base.Test
using WiltonInts84
using FixedSizeArrays

WI = WiltonInts84
p1 = Vec(1.0, 0.0, 0.0)
p2 = Vec(0.0, 1.0, 0.0)
p3 = Vec(0.0, 0.0, 0.0)
z  = Vec(0.0,0.0,1.0)

# code 01 and 00
c = p3
r, R = 0.1, 1.1
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 3
@test length(q.arcs) == 1
@test length(q.circles) == 0


# code 02
c = (p1+p2+p3)/3
r, R = 0.35, 2.0
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 6
@test length(q.arcs) == 3
@test length(q.circles) == 0

# code 10
c = (p1+p2+p3)/3
r, R = 0.1, 0.45
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 3
@test length(q.arcs) == 3
@test length(q.circles) == 1

# code 11
c = p3
r, R = 0.1, 0.3
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 2
@test length(q.arcs) == 2
@test length(q.circles) == 0

# code 12
c = Vec(0.25,0.1,0.0)
r, R = 0.2, 0.3
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 3
@test length(q.arcs) == 2
@test length(q.circles) == 0

# code 20
c = p3
r, R = 0.1, 0.9
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 3
@test length(q.arcs) == 3
@test length(q.circles) == 0
c = (p1+p2+p3)/3
r, R = 0.1, 0.45
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 3
@test length(q.arcs) == 3
@test length(q.circles) == 1

# code 21 is impossible

# code 22
c = (p1+p2+p3)/3
r, R = 0.35, 0.45
q = WI.contour(p1,p2,p3,c,r,R)
@test length(q.segments) == 6
@test length(q.arcs) == 6
@test length(q.circles) == 0
