function plot_domain(
    domain::RegularPolygon;
    num_boundary_points = 0,
    num_boundary_points_symmetry = 8,
    num_interior_points = 100,
    num_grid_points = 0,
)
    fig = GLMakie.Figure()

    axis = GLMakie.Axis(fig[1, 1], aspect = GLMakie.DataAspect())

    GLMakie.poly!(axis, vertices(domain), color = :transparent, strokewidth = 2.0)

    if num_boundary_points > 0
        GLMakie.scatter!(
            axis,
            boundary_points(domain, num_boundary_points),
            markersize = ifelse(num_boundary_points > 50, 15, 20),
        )
    end

    if num_boundary_points_symmetry > 0
        GLMakie.scatter!(
            axis,
            boundary_points_symmetry(domain, num_boundary_points_symmetry),
            markersize = ifelse(num_boundary_points_symmetry > 50, 15, 20),
        )
    end

    if num_interior_points > 0
        GLMakie.scatter!(
            axis,
            interior_points_random(domain, num_interior_points),
            markersize = ifelse(num_interior_points > 50, 15, 20),
        )
    end

    if num_grid_points > 0
        GLMakie.scatter!(
            axis,
            interior_points_grid(domain, num_grid_points)[1],
            markersize = ifelse(num_grid_points > 50, 15, 20),
        )
    end

    return fig
end

function plot_eigenfunction(u::Eigenfunction, λ; num_grid_points = 32)
    fig = GLMakie.Figure()

    axis = GLMakie.Axis(fig[1, 1], aspect = GLMakie.DataAspect())

    GLMakie.hidespines!(axis)
    GLMakie.hidedecorations!(axis)

    GLMakie.poly!(axis, vertices(u.domain), color = :transparent, strokewidth = 5.0)

    points, inside = interior_points_grid(u.domain, num_grid_points)

    values = OhMyThreads.tmap((xy, in) -> in ? u(xy, λ) : missing, points, inside)

    colorrange = let v = abs(u(center(u.domain), λ))
        (-v, v)
    end

    GLMakie.heatmap!(
        axis,
        getindex.(points, 1),
        getindex.(points, 2),
        values,
        interpolate = true;
        colorrange,
    )

    return fig
end

function plot_boundary(
    u::Eigenfunction,
    λ;
    num_boundary_points = 128,
    all_boundaries = false,
    full_boundary = false,
)
    fig = GLMakie.Figure()

    axis = GLMakie.Axis(fig[1, 1])

    if full_boundary
        points = boundary_points(u.domain, num_boundary_points)

        values = OhMyThreads.tmap(xy -> u(xy, λ), points)

        GLMakie.lines!(axis, values)

        return fig
    end

    points = boundary_points_symmetry(u.domain, num_boundary_points)

    values = OhMyThreads.tmap(xy -> u(xy, λ), points)

    GLMakie.scatter!(axis, getindex.(points, 2), values)

    points = boundary_points(u.domain, u.domain.N, num_boundary_points)

    values = OhMyThreads.tmap(xy -> u(xy, λ), points)

    GLMakie.lines!(axis, getindex.(points, 2), values)

    if all_boundaries
        for i = 1:u.domain.N-1
            GLMakie.lines!(
                axis,
                getindex.(points, 2),
                OhMyThreads.tmap(
                    xy -> u(xy, λ),
                    boundary_points(u.domain, i, num_boundary_points),
                ),
            )
        end
    end

    return fig
end

function plot_convergence!(u::Eigenfunction, λ, Ns;)
    fig = GLMakie.Figure()

    axis = GLMakie.Axis(fig[1, 1], yscale = log10)

    defects = map(Ns) do N
        sigma!(u, λ, N)
        maximum_boundary_estimate(u, λ) / norm_lower_estimate(u, λ)
    end

    GLMakie.scatterlines!(Ns, defects)

    return fig
end

function plot_convergence_vertices(λs, Ns, n)
    fig = GLMakie.Figure()

    axis = GLMakie.Axis(fig[1, 1])

    defects = map(λs, Ns) do λ, N
        u = Eigenfunction(RegularPolygon{eltype(λ)}(N))
        sigma!(u, λ, n)
        maximum_boundary_estimate(u, λ) / norm_lower_estimate(u, λ)
    end

    GLMakie.scatterlines!(Ns, defects)

    return fig
end
