struct Polar{T}
    r::T
    φ::T
end

function Polar(xy::Point2{T}) where {T}
    r = hypot(xy[1], xy[2])
    φ = atan(xy[2], xy[1])
    return Polar{T}(r, φ)
end

function Polar(xy::Point2{T}, rotation::T) where {T}
    r = hypot(xy[1], xy[2])
    φ = atan(xy[2], xy[1]) + rotation
    return Polar{T}(r, φ)
end
