type IntegrationPath{T,P}
  center::P #
  inner_radius::T
  outer_radius::T
  segments::Vector{Tuple{P,P}}
  arcs::Vector{Tuple{P,P,Int}}
  circles::Vector{Int}
  normal::P
  height::T #
  projected_center::P
  plane_inner_radius::T
  plane_outer_radius::T
end

const CLOCKWISE = -1
const COUNTERCLOCKWISE = +1


function contour(p1, p2, p3, center, inner_radius, outer_radius)

  deci = eltype(p1)
  Vertex = typeof(p1)

  t_inner = zeros(deci,2)
  t_outer = zeros(deci,2)

  inner_in  = Vector{Vertex}(3)
  inner_out = Vector{Vertex}(3)
  outer_in  = Vector{Vertex}(3)
  outer_out = Vector{Vertex}(3)

  inner_in_edge  = zeros(Int,3)
  inner_out_edge = zeros(Int,3)
  outer_in_edge  = zeros(Int,3)
  outer_out_edge = zeros(Int,3)

  is_not_empty = number_of_segments = number_of_arcs = number_of_circles = 0
  inner = outer = inner_inc = inner_outc = outer_inc = outer_outc = 0

  # centra en stralen berekenen
  inner_radius = (inner_radius > 0.0) ? inner_radius : 0.0
  outer_radius = (outer_radius > 0.0) ? outer_radius : 0.0

  # Vertex normal = f.normal()
  normal = (p1-p3) × (p2-p3)
  normal /= norm(normal)

  d = (center - p1) ⋅ normal
  absd = abs(d)
  plane_center = center - dot(normal, center - p1) * normal

  discr = inner_radius*inner_radius - d * d
  plane_inner_radius = (discr > 0.0) ? sqrt(discr) : -1.0
  discr = outer_radius*outer_radius - d * d
  plane_outer_radius = (discr > 0.0) ? sqrt(discr) : -1.0

  segments = Vector{Tuple{Vertex,Vertex}}()
  arcs = Vector{Tuple{Vertex,Vertex,Int}}()
  circles = Vector{Int}()

  # loop over de edges om segmenten te vinden
  first = zero(Vertex)
  last  = zero(Vertex)
  F = (p1,p2,p3)
  L = (p2,p3,p1)
  for i in 1:3
    first, last = F[i], L[i]
    z_inner, t_inner = incidenceLineWithSphere(center, first, last, inner_radius)
    z_outer, t_outer = incidenceLineWithSphere(center, first, last, outer_radius)
    code = 10 * z_outer + z_inner

    if code == 00
      d1 = norm(first - center)
      d2 = norm(last - center)
      D = 0.5 * (d1 + d2)
      if !(D < inner_radius || D > outer_radius)
        push!(segments, (first,last))
      end

    elseif code == 01
      t1 = t_inner[1]
      v1 = (1-t1) * first + t1 * last
      d1 = norm(first - center)
      d2 = norm(last - center)
      if d1 > d2
        push!(segments, (first, v1))
        inner_inc += 1
        inner_in[inner_inc] = v1
        inner_in_edge[inner_inc] = i
      else
        push!(segments, (v1, last))
        inner_outc += 1
        inner_out[inner_outc] = v1
        inner_out_edge[inner_outc] = i
      end

    elseif code == 02
      t1 = t_inner[1]
      t2 = t_inner[2]
      v1 = (1-t1) * first +  t1 * last
      inner_inc += 1
      inner_in[inner_inc] = v1
      inner_in_edge[inner_inc] = i
      v2 = (1-t2) * first + t2 * last
      inner_outc += 1
      inner_out[inner_outc] = v2
      inner_out_edge[inner_outc] = i
      push!(segments, (first,v1))
      push!(segments, (v2,last))

    elseif code == 10
      t1 = t_outer[1]
      v1 = (1-t1) * first +  t1 * last
      d1 = norm(first - center)
      d2 = norm(last - center)
      if d1 < d2
          push!(segments, (first, v1))
          outer_inc += 1
          outer_in[outer_inc] = v1
          outer_in_edge[outer_inc] = i
      else
          push!(segments, (v1, last))
          outer_outc += 1
          outer_out[outer_outc] = v1
          outer_out_edge[outer_outc] = i
      end


    elseif code == 11
      t1 = t_inner[1]
      t2 = t_outer[1]
      v1 = (1-t1) * first + t1 * last
      v2 = (1-t2) * first + t2 * last
      if t1 < t2
        push!(segments, (v1, v2))
        inner_outc += 1
        inner_out[inner_outc] = v1
        inner_out_edge[inner_outc] = i
        outer_inc += 1
        outer_in[outer_inc] = v2
        outer_in_edge[outer_inc] = i
      else
        push!(segments, (v2, v1))
        inner_inc += 1
        inner_in[inner_inc] = v1
        inner_in_edge[inner_inc] = i
        outer_outc += 1
        outer_out[outer_outc] = v2
        outer_out_edge[outer_outc] = i
      end

    elseif code == 12
      t1 = t_inner[1]
      t2 = t_inner[2]
      t3 = t_outer[1]
      v1 = (1-t1) * first + t1 * last
      v2 = (1-t2) * first + t2 * last
      v3 = (1-t3) * first + t3 * last
      if t1 < t3
        push!(segments, (first, v1))
        push!(segments, (v2, v3))
        inner_inc += 1
        inner_in[inner_inc] = v1
        inner_in_edge[inner_inc] = i
        inner_outc += 1
        inner_out[inner_outc] = v2
        inner_out_edge[inner_outc] = i
        outer_inc += 1
        outer_in[outer_inc] = v3
        outer_in_edge[outer_inc] = i
      else
        push!(segments, (v3, v1))
        push!(segments, (v2, last))
        inner_inc += 1
        inner_in[inner_inc] = v1
        inner_in_edge[inner_inc] = i
        inner_outc += 1
        inner_out[inner_outc] = v2
        inner_out_edge[inner_outc] = i
        outer_outc += 1
        outer_out[outer_outc] = v3
        outer_out_edge[outer_outc] = i
      end

    elseif code == 20
      t1 = t_outer[1]
      t2 = t_outer[2]
      v1 = (1-t1) * first + t1 * last
      outer_outc += 1
      outer_out[outer_outc] = v1
      outer_out_edge[outer_outc] = i
      v2 = (1-t2) * first + t2 * last
      outer_inc += 1
      outer_in[outer_inc] = v2
      outer_in_edge[outer_inc] = i
      push!(segments, (v1, v2))

    elseif code == 22
      t1 = t_inner[1]
      t2 = t_inner[2]
      t3 = t_outer[1]
      t4 = t_outer[2]
      v1 = (1-t1) * first + t1 * last
      v2 = (1-t2) * first + t2 * last
      v3 = (1-t3) * first + t3 * last
      v4 = (1-t4) * first + t4 * last
      inner_inc += 1
      inner_in[inner_inc] = v1
      inner_in_edge[inner_inc] = i
      inner_outc += 1
      inner_out[inner_outc] = v2
      inner_out_edge[inner_outc] = i
      outer_outc += 1
      outer_out[outer_outc] = v3
      outer_out_edge[outer_outc] = i
      outer_inc += 1
      outer_in[outer_inc] = v4
      outer_in_edge[outer_inc] = i
      push!(segments, (v3, v1))
      push!(segments, (v2, v4))
    end

  end # end for i

  # bereken het aantal snijpunten met de twee cirkels
  inner = inner_inc + inner_outc
  outer = outer_inc + outer_outc

  @assert inner_inc == inner_outc
  @assert outer_inc == outer_outc

  # construct the inner arcs
  if inner == 0
  elseif inner == 2
    push!(arcs, (inner_in[1], inner_out[1], CLOCKWISE))
  elseif inner == 4
    if inner_in_edge[1] == inner_out_edge[1]
      push!(arcs, (inner_in[1],inner_out[2],CLOCKWISE))
      push!(arcs, (inner_in[2], inner_out[1], CLOCKWISE))
    elseif inner_in_edge[1] == inner_out_edge[2]
      push!(arcs, (inner_in[1],inner_out[1],CLOCKWISE))
      push!(arcs, (inner_in[2], inner_out[2],CLOCKWISE))
    else
      ed1 = inner_in_edge[1]
      ed2 = mod1(ed1-1, 3)
      if inner_out_edge[1] == ed2
        # This branch -I think- cannot be reached
        push!(arcs, (inner_in[1], inner_out[1], CLOCKWISE))
        push!(arcs, (inner_in[2], inner_out[2], CLOCKWISE))
      else
        push!(arcs, (inner_in[1], inner_out[2], CLOCKWISE))
        push!(arcs, (inner_in[2], inner_out[1], CLOCKWISE))
      end
    end
  elseif inner == 6
    for i in 1:3
      ed1 = inner_in_edge[i]
      ed2 = mod1(ed1-1,3)
      for j in 1:3
        if ed2 == inner_out_edge[j]
          push!(arcs, (inner_in[i], inner_out[j], CLOCKWISE))
        end
      end
    end
  end

  # construct the outer arcs
  if outer == 0
  elseif outer == 2
    push!(arcs, (outer_in[1], outer_out[1], COUNTERCLOCKWISE))
  elseif outer == 4
    if outer_in_edge[1] == outer_out_edge[1]
      push!(arcs, (outer_in[1], outer_out[2], COUNTERCLOCKWISE))
      push!(arcs, (outer_in[2], outer_out[1], COUNTERCLOCKWISE))
    elseif outer_in_edge[1] == outer_out_edge[2]
      # This branch cannot be reached
      push!(arcs, (outer_in[1], outer_out[1], COUNTERCLOCKWISE))
      push!(arcs, (outer_in[2], outer_out[2], COUNTERCLOCKWISE))
    else
      ed1 = outer_in_edge[1]
      ed2 = mod1(ed1+1,3)
      if outer_out_edge[1] == ed2
        push!(arcs, (outer_in[1], outer_out[1], COUNTERCLOCKWISE))
        push!(arcs, (outer_in[2], outer_out[2], COUNTERCLOCKWISE))
      else
        push!(arcs, (outer_in[1], outer_out[2], COUNTERCLOCKWISE))
        push!(arcs, (outer_in[2], outer_out[1], COUNTERCLOCKWISE))
      end
    end
  elseif outer == 6
    for i in 1:3
      ed1 = outer_in_edge[i]
      ed2 = mod1(ed1+1, 3)
      for j in 1:3
        if ed2==outer_out_edge[j]
          push!(arcs, (outer_in[i], outer_out[j], COUNTERCLOCKWISE))
        end
      end
    end
  end

  # construct the circle contributions
  if inside(plane_center,p1,p2,p3,normal)
    d1 = distancetoline(plane_center, p1, p2)
    d2 = distancetoline(plane_center, p2, p3)
    d3 = distancetoline(plane_center, p3, p1)
    if d1>plane_inner_radius && d2>plane_inner_radius && d3>plane_inner_radius && plane_inner_radius>0
      push!(circles, CLOCKWISE)
    end
    if d1>plane_outer_radius && d2>plane_outer_radius && d3>plane_outer_radius && plane_outer_radius>0
      push!(circles, COUNTERCLOCKWISE)
    end
  end

  IntegrationPath(
    center,
    inner_radius,
    outer_radius,
    segments,
    arcs,
    circles,
    normal,
    d,
    plane_center,
    plane_inner_radius,
    plane_outer_radius
  )
end


function inside(v,p1,p2,p3,n)
  dot((p2-p1)×(v-p1), n) <= 0 && return false
  dot((p3-p2)×(v-p2), n) <= 0 && return false
  dot((p1-p3)×(v-p3), n) <= 0 && return false
  return true
end


function distancetoline(p,a,b)
  t = b-a
  t /= norm(t)
  u = p-a
  norm(u - dot(u,t)*t)
end


"""
  r, t = incidenceLineWithSphere(v, first, last, r)

Returns the number of crossings and their barycentric coordinates for a circle
with center v and radius r and a segment [a,b]. If the circle is tangent to the
segment, zero crossings are reported.
"""
function incidenceLineWithSphere(v, first, last, r)

  result = 0
  t = fill(zero(eltype(v)), 2)

  r<=0 && return result,t

  a = dot((last - first) , (last - first))
  b =  2 * dot( (first - v) , (last - first))
  c = dot((first - v) , (first - v)) - r*r
  d = b * b - 4 * a * c
  d < 0 && return result,t

  f = sqrt(d)
  t[result+1] = (-b - f) / (2 * a)
  t[result+1]>0 && t[result+1]<=1 && (result+=1)
  t[result+1] = (-b + f) / (2 * a)
  t[result+1]>0 && t[result+1]<=1 && (result+=1)
  return result,t

end
