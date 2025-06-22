# Set environment variable to avoid graphics issues in some environments
ENV["GKSwstype"]="100"

println("Loading required packages and modules...")
include("../src/CollocationPoints.jl")
using .CollocationPoints 

using CairoMakie
using LinearAlgebra

# Create figures directory if it doesn't exist
const FIG_DIR = "figures"
if !isdir(FIG_DIR)
    mkdir(FIG_DIR)
end

# Define a function to create a clean axis without decorations
function create_clean_axis(fig)
    ax = Axis(fig[1, 1], aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    return ax
end

function save_figure(name, fig)
    path = joinpath(FIG_DIR, name)
    save(path, fig, px_per_unit = 2)
    println("Plot saved to $path")
end

println("Setup complete. Running all examples...")

# --- Example 1: Uniform points in a circle ---
println("\nGenerating uniform points in a circle...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = zeros(0, 2)
    local radius_fn(p, ctps) = fill(0.05, size(p, 1))
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 4, strokewidth = 0)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :blue, linewidth = 2)
    save_figure("uniform_points_in_circle.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for uniform circle: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 2: Variable density points in a circle (center focus) ---
println("\nGenerating variable density points in a circle (center focus)...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = [0.0 0.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.08 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 4, strokewidth = 0)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :blue, linewidth = 2)
    scatter!(ax, ctps[:, 1], ctps[:, 2], color = :red, markersize = 8, marker = :circle)
    save_figure("variable_density_circle_center.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for variable density circle (center): $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 3: Variable density points in a circle (line focus) ---
println("\nGenerating variable density points in a circle (line focus)...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = hcat(range(-0.5, 0.5, length=10), zeros(10))
    local radius_fn(p, ctps) = 0.005 .+ 0.08 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 4, strokewidth = 0)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :blue, linewidth = 2)
    scatter!(ax, ctps[:, 1], ctps[:, 2], color = :red, markersize = 6)
    lines!(ax, ctps[:, 1], ctps[:, 2], color = :red, linewidth = 2)
    save_figure("variable_density_circle_line.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for variable density circle (line): $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 4: Variable density points in a circle (two-point focus) ---
println("\nGenerating variable density points in a circle (two-point focus)...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_circ(0.0, 0.0, 1.0, 2 / hbdy)
    local ctps = [-0.75 0.0; 0.75 0.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.08 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 4, strokewidth = 0)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :blue, linewidth = 2)
    scatter!(ax, ctps[:, 1], ctps[:, 2], color = :red, markersize = 8, marker = :circle)
    save_figure("variable_density_circle_two_points.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for variable density circle (two points): $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 5: Uniform points in a rectangle ---
println("\nGenerating uniform points in a rectangle...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.05
    local ptol = 0.01
    local bdy, fd = draw_rect(-1.0, -1.0, 1.0, 1.0, 2 / hbdy)
    local ctps = zeros(0, 2)
    local radius_fn(p, ctps) = fill(0.05, size(p, 1))
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 4, strokewidth = 0)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :blue, linewidth = 2)
    save_figure("uniform_points_in_rectangle.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for uniform rectangle: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 6: Variable density points in a rectangle ---
println("\nGenerating variable density points in a rectangle...")
try
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.02
    local ptol = 0.001
    local bdy, fd = draw_rect(-1.0, -1.0, 1.0, 1.0, 2 / hbdy)
    local ctps = [-0.5 1.0; 0.5 1.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.05 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 4, strokewidth = 0)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :blue, linewidth = 2)
    scatter!(ax, ctps[:, 1], ctps[:, 2], color = :red, markersize = 8, marker = :circle)
    save_figure("variable_density_rectangle.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for variable density rectangle: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 7: L-shaped domain with corner singularity ---
println("\nGenerating points in an L-shaped domain with corner singularity...")
try
    local demo_file = joinpath("demos", "Lshape.txt")
    local box = [-1.0, 1.0, -1.0, 1.0]
    local hbdy = 0.025
    local ptol = 0.001
    local b_xy, fd = make_domain(demo_file)
    local bdy = bsmooth(b_xy, hbdy)
    local ctps = [0.0 0.0]
    local radius_fn(p, ctps) = 0.005 .+ 0.05 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 4, strokewidth = 0, opacity = 0.8)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :blue, linewidth = 2)
    scatter!(ax, ctps[:, 1], ctps[:, 2], color = :red, markersize = 8, marker = :circle)
    save_figure("l_shape_domain_with_singularity.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for L-shaped domain: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 8: Lake model with boundary density ---
println("\nGenerating points in a lake model with boundary density...")
try
    local demo_file = joinpath("demos", "lake.txt")
    local box = [100.0, 634.0, 145.0, 799.0]
    local hbdy = 5.0
    local ptol = 1.0
    local b_xy, fd = make_domain(demo_file)
    local bdy = bsmooth(b_xy, hbdy)
    local ctps = bdy
    local radius_fn(p, ctps) = 2.0 .+ 0.1 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 3.5, strokewidth = 0)
    lines!(ax, bdy[:, 1], bdy[:, 2], color = :navy, linewidth = 2.5)
    save_figure("lake_model_with_boundary_density.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for lake model: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

# --- Example 9: Island with an inner lake ---
println("\nGenerating points for an island with an inner lake...")
try
    local lake_file = joinpath("demos", "lake.txt")
    local island_file = joinpath("demos", "island.txt")
    local box = [100.0, 634.0, 145.0, 799.0]
    local hbdy = 3.0
    local ptol = 1.0
    local b1_xy, fd1 = make_domain(lake_file)
    local b2_xy, fd2 = make_domain(island_file)
    local fd(p) = max.(fd1(p), -fd2(p))
    local b1_bdy = bsmooth(b1_xy, hbdy)
    local b2_bdy = bsmooth(b2_xy, hbdy)
    local bdy = vcat(b1_bdy, b2_bdy)
    local ctps = bdy
    local radius_fn(p, ctps) = 2.0 .+ 0.2 .* min_pdist2(p, ctps)
    local xy = CollocationPoints2D(fd, box, ctps, ptol, radius_fn)
    
    fig = Figure(size = (800, 800))
    ax = create_clean_axis(fig)
    scatter!(ax, xy[:, 1], xy[:, 2], color = :black, markersize = 3.5, strokewidth = 0)
    lines!(ax, b1_bdy[:, 1], b1_bdy[:, 2], color = :navy, linewidth = 2.5)
    lines!(ax, b2_bdy[:, 1], b2_bdy[:, 2], color = :darkgreen, linewidth = 2.5)
    save_figure("island_with_inner_lake.png", fig)
    println("Generated $(size(xy, 1)) nodes.")
catch e
    println("Error generating points for island model: $e")
    showerror(stdout, e, catch_backtrace())
    println()
end

println("\nAll examples finished.")