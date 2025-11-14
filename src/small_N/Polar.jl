"""
    Polar{T}(r::T, φ::T)
    Polar(r::T, φ::T)
    Polar(xy::Point2)
    Polar(xy::Point2, rotation)

Represents a point in polar coordinates.

The constructor taking `rotation` as an argument gives the polar
coordinates of the point `xy` rotated by `rotation` in the positive
direction.
"""
struct Polar{T}
    r::T
    φ::T
end

Polar(r::T, φ::T) where {T<:Number} = Polar{T}(r, φ)

function Polar(xy::Point2)
    r = hypot(xy[1], xy[2])
    φ = atan(xy[2], xy[1])
    return Polar(r, φ)
end

function Polar(xy::Point2, rotation)
    # Rotate input before conversion to polar coordinates. In
    # principle it should be faster to rotate after conversion, but it
    # is more difficult to ensure that φ is in the interval [-π, π].
    s, c = sincos(rotation)
    x = c * xy[1] - s * xy[2]
    y = s * xy[1] + c * xy[2]

    r = hypot(x, y)
    φ = atan(y, x)

    return Polar(r, φ)
end
