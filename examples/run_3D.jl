# example/run_3D.jl
ENV["GKSwstype"]="100"

include("../src/CollocationPoints3.jl")
using .CollocationPoints3

using CairoMakie
using LinearAlgebra

const FIG_DIR = "figures_3d"
if !isdir(FIG_DIR)
    mkdir(FIG_DIR)
end

function save_figure(name, fig)
    path = joinpath(FIG_DIR, name)
    save(path, fig)
    println("Plot saved to $path")
end

println("Generating uniform points in a sphere...")

bbox = [-0.5, 0.5, -0.5, 0.5, -0.5, 0.5]
radius = 0.5
ptol = 0.01

_, sdf_sphere = draw_sphere(0.0, 0.0, 0.0, radius, 10)

ctps = zeros(0, 3)

# Use a moderate spacing for good results
spacing = 0.04
radius_fn(p, ctps) = radius_3d_unif(p, spacing)

println("Starting sphere generation...")
xyz = CollocationPoints3D(sdf_sphere, bbox, ctps, ptol, radius_fn)
println("Generated $(size(xyz, 1)) nodes in sphere.")

fig = Figure(size = (800, 800))
ax = Axis3(fig[1, 1], aspect = :data, 
           xlabel = "x", ylabel = "y", zlabel = "z")

scatter!(ax, xyz[:, 1], xyz[:, 2], xyz[:, 3], 
         color = :black, markersize = 5, label = " ")



axislegend(ax)

save_figure("fixed_3d_point_cloud.png", fig)
