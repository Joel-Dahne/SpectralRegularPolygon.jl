struct Polar{T}
    r::T
    φ::T
end

Polar(r::T, φ::T) where {T<:Number} = Polar{T}(r, φ)

function Polar(xy::Point2)
    r = sqrt(xy[1]^2 + xy[2]^2)
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

    r = sqrt(x^2 + y^2)
    φ = atan(y, x)

    return Polar(r, φ)
end
