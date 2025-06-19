# Set environment variable to avoid graphics issues
ENV["GKSwstype"]="100"

println("Loading required packages and modules...")
using CairoMakie
include("CollocationPoints.jl")
using .CollocationPoints
using LinearAlgebra

# Create figures directory if it doesn't exist
if !isdir("figures")
    mkdir("figures")
end

println("Setup complete. Running all examples...")

# Define a function to create a clean axis without decorations
function create_clean_axis(fig)
    ax = Axis(fig[1, 1], aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

println("\nGenerating uniform points in a circle...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = zeros(0, 2)
    local radius_fn(p, ctps) = fill(0.05, size(p, 1))
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 4, 
        strokewidth = 0
    )
    
    # Plot boundary
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :blue, 
        linewidth = 2
    )
    
    # Save figure
    save("figures/uniform_points_in_circle.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/uniform_points_in_circle.png")
catch e
    println("Error generating points for uniform circle: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating variable density points in a circle (center focus)...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = [0.0 0.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.08 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 4, 
        strokewidth = 0
    )
    
    # Plot boundary
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :blue, 
        linewidth = 2
    )
    
    # Plot the focus point
    scatter!(ax, ctps[:, 1], ctps[:, 2], 
        color = :red, 
        markersize = 8, 
        marker = :circle
    )
    
    # Save figure
    save("figures/variable_density_circle_center.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/variable_density_circle_center.png")
catch e
    println("Error generating points for variable density circle (center): $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating variable density points in a circle (line focus)...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = hcat(range(-0.5, 0.5, length=10), zeros(10))
    local radius_fn(p, ctps) = 0.005 .+ 0.08 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 4, 
        strokewidth = 0
    )
    
    # Plot boundary
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :blue, 
        linewidth = 2
    )
    
    # Plot the focus points
    scatter!(ax, ctps[:, 1], ctps[:, 2], 
        color = :red, 
        markersize = 6
    )
    
    # Connect focus points with a line
    lines!(ax, ctps[:, 1], ctps[:, 2], 
        color = :red, 
        linewidth = 2
    )
    
    # Save figure
    save("figures/variable_density_circle_line.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/variable_density_circle_line.png")
catch e
    println("Error generating points for variable density circle (line): $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating variable density points in a circle (two-point focus)...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = [-0.75 0.0; 0.75 0.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.08 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 4, 
        strokewidth = 0
    )
    
    # Plot boundary
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :blue, 
        linewidth = 2
    )
    
    # Plot the focus points
    scatter!(ax, ctps[:, 1], ctps[:, 2], 
        color = :red, 
        markersize = 8, 
        marker = :circle
    )
    
    # Save figure
    save("figures/variable_density_circle_two_points.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/variable_density_circle_two_points.png")
catch e
    println("Error generating points for variable density circle (two points): $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating uniform points in a rectangle...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.05
    local ptol = 0.01
    local bdy, fd = draw_rect(-1.0, -1.0, 1.0, 1.0, 2 / hbdy)
    local ctps = zeros(0, 2)
    local radius_fn(p, ctps) = fill(0.05, size(p, 1))
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 4, 
        strokewidth = 0
    )
    
    # Plot boundary
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :blue, 
        linewidth = 2
    )
    
    # Save figure
    save("figures/uniform_points_in_rectangle.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/uniform_points_in_rectangle.png")
catch e
    println("Error generating points for uniform rectangle: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating variable density points in a rectangle...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_rect(-1.0, -1.0, 1.0, 1.0, 2 / hbdy)
    local ctps = [-0.5 1.0; 0.5 1.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.05 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 4, 
        strokewidth = 0
    )
    
    # Plot boundary
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :blue, 
        linewidth = 2
    )
    
    # Plot the focus points
    scatter!(ax, ctps[:, 1], ctps[:, 2], 
        color = :red, 
        markersize = 8, 
        marker = :circle
    )
    
    # Save figure
    save("figures/variable_density_rectangle.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/variable_density_rectangle.png")
catch e
    println("Error generating points for variable density rectangle: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating points in an L-shaped domain with corner singularity...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.025
    local ptol = 0.001
    local b_xy, fd = make_domain("demos/Lshape.txt")
    local bdy = bsmooth(b_xy, hbdy)
    local ctps = [0.0 0.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.05 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 4, 
        strokewidth = 0,
        opacity = 0.8
    )
    
    # Plot boundary as a smooth line
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :blue, 
        linewidth = 2
    )
    
    # Mark the singularity point
    scatter!(ax, ctps[:, 1], ctps[:, 2], 
        color = :red, 
        markersize = 8, 
        marker = :circle
    )
    
    # Save figure
    save("figures/l_shape_domain_with_singularity.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/l_shape_domain_with_singularity.png")
catch e
    println("Error generating points for L-shaped domain: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating points in a lake model with boundary density...")
try
    local box = [100.0, 634.0, 145.0, 799.0]
    local hbdy = 5.0
    local ptol = 1.0
    local b_xy, fd = make_domain("demos/lake.txt")
    local bdy = bsmooth(b_xy, hbdy)
    local ctps = bdy
    local radius_fn(p, ctps) = 2.0 .+ 0.1 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points only - no polygons or lines inside
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 3.5, 
        strokewidth = 0
    )
    
    # Plot only the boundary line
    lines!(ax, bdy[:, 1], bdy[:, 2], 
        color = :navy, 
        linewidth = 2.5
    )
    
    # Save figure
    save("figures/lake_model_with_boundary_density.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/lake_model_with_boundary_density.png")
catch e
    println("Error generating points for lake model: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nGenerating points for an island with an inner lake...")
try
    local box = [100.0, 634.0, 145.0, 799.0]
    local hbdy = 3.0
    local ptol = 1.0
    local b1_xy, fd1 = make_domain("demos/lake.txt")
    local b2_xy, fd2 = make_domain("demos/island.txt")
    local fd(p) = max.(fd1(p), -fd2(p))
    local b1_bdy = bsmooth(b1_xy, hbdy)
    local b2_bdy = bsmooth(b2_xy, hbdy)
    local bdy = vcat(b1_bdy, b2_bdy)
    local ctps = bdy
    local radius_fn(p, ctps) = 2.0 .+ 0.2 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    # Create figure with Makie
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    
    # Plot points only - no polygons or interior lines
    scatter!(ax, xy[:, 1], xy[:, 2], 
        color = :black, 
        markersize = 3.5,  # Larger point size
        strokewidth = 0
    )
    
    # Plot lake boundary (outer boundary)
    lines!(ax, b1_bdy[:, 1], b1_bdy[:, 2], 
        color = :navy, 
        linewidth = 2.5
    )
    
    # Plot island boundary (inner boundary)
    lines!(ax, b2_bdy[:, 1], b2_bdy[:, 2], 
        color = :darkgreen, 
        linewidth = 2.5
    )
    
    # Save figure
    save("figures/island_with_inner_lake.png", fig, px_per_unit = 2)
    println("Generated $(size(xy, 1)) nodes. Plot saved to figures/island_with_inner_lake.png")
catch e
    println("Error generating points for island model: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nAll examples finished.")