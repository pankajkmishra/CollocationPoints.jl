module CollocationPoints3

using LinearAlgebra
using Random

export CollocationPoints3D, draw_sphere, draw_box, draw_cylinder, min_pdist3, radius_3d_unif

"""
    CollocationPoints3D(fd, bbox, ctps, ptol, radius_fn)

Generate points in a 3D domain using a fast block-filling algorithm.

# Arguments
- `fd`: A signed distance function that takes an Nx3 array of points and returns
        an N-length vector of negative values inside the domain and positive outside.
- `bbox`: The bounding box as [xmin, xmax, ymin, ymax, zmin, zmax].
- `ctps`: Control points (Nx3) to be included in the output and used for density control.
- `ptol`: Tolerance for boundary points.
- `radius_fn`: A function that takes an Nx3 array of points and returns an N-length vector
               of desired radii at those points.

# Returns
- An Nx3 array of points in the domain.
"""
function CollocationPoints3D(fd, bbox, ctps, ptol, radius_fn)
    # Generate points using improved block sampling
    points = block_sampling(fd, bbox, radius_fn, ctps)
    
    if isempty(points) 
        return zeros(0, 3) 
    end
    
    # Keep only points inside the domain
    inside_idx = fd(points) .< eps()
    points = points[inside_idx, :]
    
    # Remove points that are already in the control points
    if !isempty(ctps) && !isempty(points)
        p_tuples = [tuple(row...) for row in eachrow(points)]
        ctps_tuples = [tuple(row...) for row in eachrow(ctps)]
        p_tuples = setdiff(p_tuples, ctps_tuples)
        points = isempty(p_tuples) ? zeros(0, 3) : vcat([collect(t)' for t in p_tuples]...)
    end
    
    # Add control points
    pfix = isempty(ctps) ? zeros(0, 3) : unique(ctps, dims=1)
    points = vcat(pfix, points)
    
    if isempty(points) 
        return zeros(0, 3) 
    end
    
    # Remove boundary points
    boundary_idx = abs.(fd(points)) .< ptol
    keep_idx = .!boundary_idx
    points = points[keep_idx, :]
    
    return points
end

"""
    block_sampling(fd, bbox, radius_fn, ctps)

Generate points using block-based sampling.

# Arguments
- `fd`: Signed distance function.
- `bbox`: Bounding box [xmin, xmax, ymin, ymax, zmin, zmax].
- `radius_fn`: Function that computes the radius at each point.
- `ctps`: Control points.

# Returns
- Array of 3D points.
"""
function block_sampling(fd, bbox, radius_fn, ctps)
    xmin, xmax, ymin, ymax, zmin, zmax = bbox
    
    # Determine minimum radius by sampling
    sample_pts = [
        [xmin, ymin, zmin];
        [xmax, ymin, zmin];
        [xmin, ymax, zmin];
        [xmax, ymax, zmin];
        [xmin, ymin, zmax];
        [xmax, ymin, zmax];
        [xmin, ymax, zmax];
        [xmax, ymax, zmax];
        [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2];
    ]
    
    radii = radius_fn(sample_pts, ctps)
    min_radius = minimum(radii)
    
    # Calculate domain dimensions
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    
    # Calculate number of cells in each direction
    nx = max(2, ceil(Int, dx / min_radius))
    ny = max(2, ceil(Int, dy / min_radius))
    nz = max(2, ceil(Int, dz / min_radius))
    
    # Cell size
    cell_dx = dx / nx
    cell_dy = dy / ny
    cell_dz = dz / nz
    
    # Initialize points list
    points = []
    
    # Jitter factor for randomization
    jitter = 0.5
    
    # Generate candidate points at each cell center with jitter
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                # Calculate cell center
                x = xmin + (i - 0.5) * cell_dx
                y = ymin + (j - 0.5) * cell_dy
                z = zmin + (k - 0.5) * cell_dz
                
                # Add random jitter
                x_jit = x + (rand() - 0.5) * jitter * cell_dx
                y_jit = y + (rand() - 0.5) * jitter * cell_dy
                z_jit = z + (rand() - 0.5) * jitter * cell_dz
                
                # Keep points inside the bounding box
                x_jit = clamp(x_jit, xmin + 1e-10, xmax - 1e-10)
                y_jit = clamp(y_jit, ymin + 1e-10, ymax - 1e-10)
                z_jit = clamp(z_jit, zmin + 1e-10, zmax - 1e-10)
                
                push!(points, [x_jit, y_jit, z_jit])
            end
        end
    end
    
    # Convert to matrix
    points_matrix = isempty(points) ? zeros(0, 3) : hcat(points...)'
    
    # Filter to keep only points inside domain or near boundary
    if !isempty(points_matrix)
        # Evaluate SDF at all points
        sdf_values = fd(points_matrix)
        
        # Keep points inside or within tolerance of boundary
        inside_idx = sdf_values .<= 0.01
        points_matrix = points_matrix[inside_idx, :]
    end
    
    # Filter points to maintain spacing
    if !isempty(points_matrix)
        points_matrix = filter_by_distance(points_matrix, radius_fn, ctps)
    end
    
    return points_matrix
end

"""
    filter_by_distance(points, radius_fn, ctps)

Filter points so that they maintain proper spacing based on the radius function.

# Arguments
- `points`: Input points matrix.
- `radius_fn`: Function that computes the radius at each point.
- `ctps`: Control points.

# Returns
- Filtered points.
"""
function filter_by_distance(points, radius_fn, ctps)
    if isempty(points)
        return points
    end
    
    # Compute radius at each point
    radii = radius_fn(points, ctps)
    
    # Initialize filter
    keep = trues(size(points, 1))
    
    # Filter points, keeping a subset that maintains proper spacing
    for i in 1:size(points, 1)
        if !keep[i] 
            continue
        end
        
        p1 = points[i, :]
        r1 = radii[i]
        
        # Check against all other kept points
        for j in i+1:size(points, 1)
            if !keep[j] 
                continue
            end
            
            p2 = points[j, :]
            r2 = radii[j]
            
            # Combine radii for check
            r_check = (r1 + r2) / 2
            
            # Check distance
            dist = norm(p1 - p2)
            
            # If too close, discard the second point
            if dist < r_check
                keep[j] = false
            end
        end
    end
    
    # Return filtered points
    return points[keep, :]
end

"""
    min_pdist3(points, ctps)

Calculate the minimum distance from each point in 'points' to any point in 'ctps'.

# Arguments
- `points`: An Nx3 array of query points.
- `ctps`: An Mx3 array of control points.

# Returns
- An N-length vector of minimum distances.
"""
function min_pdist3(points, ctps)
    if isempty(ctps) || isempty(points)
        return fill(Inf, size(points, 1))
    end
    
    dists = zeros(size(points, 1))
    for i in 1:size(points, 1)
        p = points[i, :]
        dists[i] = minimum(sqrt.(sum((ctps .- p').^2, dims=2)))
    end
    
    return dists
end

"""
    radius_3d_unif(xyz, rFactor)

Return uniform grain radius at all locations.

# Arguments
- `xyz`: Nx3 array of points.
- `rFactor`: Constant radius value.

# Returns
- Vector of uniform radii.
"""
function radius_3d_unif(xyz, rFactor)
    return fill(rFactor, size(xyz, 1))
end

"""
    draw_sphere(xc, yc, zc, r, n)

Create a sphere with center (xc, yc, zc) and radius r.

# Arguments
- `xc`, `yc`, `zc`: Center coordinates of the sphere.
- `r`: Radius of the sphere.
- `n`: Number of points to use for approximating the sphere surface.

# Returns
- A tuple (surface, sdf) where surface is a matrix of points on the sphere surface
  and sdf is the signed distance function for the sphere.
"""
function draw_sphere(xc, yc, zc, r, n)
    n_int = round(Int, sqrt(n))
    
    theta = range(0, pi, length=n_int)
    phi = range(0, 2*pi, length=2*n_int)
    
    points = []
    for t in theta
        for p in phi
            x = r * sin(t) * cos(p) + xc
            y = r * sin(t) * sin(p) + yc
            z = r * cos(t) + zc
            push!(points, [x, y, z])
        end
    end
    
    surface = length(points) > 0 ? hcat(points...)' : zeros(0, 3)
    
    # SDF function for sphere (properly handles matrix input)
    function sdf(p)
        if size(p, 2) != 3
            error("Expected points matrix with 3 columns, got size $(size(p))")
        end
        return sqrt.((p[:,1] .- xc).^2 .+ (p[:,2] .- yc).^2 .+ (p[:,3] .- zc).^2) .- r
    end
    
    return (surface, sdf)
end

"""
    draw_box(x1, y1, z1, x2, y2, z2, n)

Create a box with corners at (x1, y1, z1) and (x2, y2, z2).

# Arguments
- `x1`, `y1`, `z1`: Minimum corner coordinates of the box.
- `x2`, `y2`, `z2`: Maximum corner coordinates of the box.
- `n`: Approximate number of points to use for representing the box surface.

# Returns
- A tuple (surface, sdf) where surface is a matrix of points on the box surface
  and sdf is the signed distance function for the box.
"""
function draw_box(x1, y1, z1, x2, y2, z2, n)
    n_side = max(10, round(Int, (n/6)^(1/2)))
    
    pts = []
    
    for x in range(x1, stop=x2, length=n_side)
        for y in range(y1, stop=y2, length=n_side)
            push!(pts, [x, y, z1])
        end
    end
    
    for x in range(x1, stop=x2, length=n_side)
        for y in range(y1, stop=y2, length=n_side)
            push!(pts, [x, y, z2])
        end
    end
    
    for x in range(x1, stop=x2, length=n_side)
        for z in range(z1, stop=z2, length=n_side)[2:end-1]
            push!(pts, [x, y1, z])
        end
    end
    
    for x in range(x1, stop=x2, length=n_side)
        for z in range(z1, stop=z2, length=n_side)[2:end-1]
            push!(pts, [x, y2, z])
        end
    end
    
    for y in range(y1, stop=y2, length=n_side)[2:end-1]
        for z in range(z1, stop=z2, length=n_side)[2:end-1]
            push!(pts, [x1, y, z])
        end
    end
    
    for y in range(y1, stop=y2, length=n_side)[2:end-1]
        for z in range(z1, stop=z2, length=n_side)[2:end-1]
            push!(pts, [x2, y, z])
        end
    end
    
    surface = length(pts) > 0 ? hcat(pts...)' : zeros(0, 3)
    
    # SDF function for box
    function sdf(p)
        if size(p, 2) != 3
            error("Expected points matrix with 3 columns, got size $(size(p))")
        end
        
        dx = maximum(hcat(x1 .- p[:, 1], p[:, 1] .- x2), dims=2)
        dy = maximum(hcat(y1 .- p[:, 2], p[:, 2] .- y2), dims=2)
        dz = maximum(hcat(z1 .- p[:, 3], p[:, 3] .- z2), dims=2)
        return maximum(hcat(dx, dy, dz), dims=2)[:]
    end
    
    return (surface, sdf)
end

"""
    draw_cylinder(x1, y1, z1, x2, y2, z2, r, n)

Create a cylinder with axis from point (x1, y1, z1) to (x2, y2, z2) and radius r.

# Arguments
- `x1`, `y1`, `z1`: First endpoint of the cylinder axis.
- `x2`, `y2`, `z2`: Second endpoint of the cylinder axis.
- `r`: Radius of the cylinder.
- `n`: Approximate number of points to use for representing the cylinder surface.

# Returns
- A tuple (surface, sdf) where surface is a matrix of points on the cylinder surface
  and sdf is the signed distance function for the cylinder.
"""
function draw_cylinder(x1, y1, z1, x2, y2, z2, r, n)
    axis = [x2 - x1, y2 - y1, z2 - z1]
    len = norm(axis)
    axis = axis / len
    
    if abs(axis[3]) < 0.9
        u = normalize(cross(axis, [0, 0, 1]))
    else
        u = normalize(cross(axis, [1, 0, 0]))
    end
    v = cross(axis, u)
    
    n_circ = max(16, round(Int, sqrt(n * r / len)))
    n_len = max(10, round(Int, sqrt(n * len / r)))
    
    theta = range(0, 2π, length=n_circ+1)[1:n_circ]
    z_vals = range(0, len, length=n_len)
    
    points = []
    for t in theta
        for z in z_vals
            p = [x1, y1, z1] + z * axis + r * (cos(t) * u + sin(t) * v)
            push!(points, p)
        end
    end
    
    for i in 1:n_circ
        t = 2π * (i-1)/n_circ
        for rad in range(0, r, length=5)[2:end]
            p1 = [x1, y1, z1] + rad * (cos(t) * u + sin(t) * v)
            push!(points, p1)
            
            p2 = [x2, y2, z2] + rad * (cos(t) * u + sin(t) * v)
            push!(points, p2)
        end
    end
    
    surface = length(points) > 0 ? hcat(points...)' : zeros(0, 3)
    
    # SDF function for cylinder
    function sdf(p)
        if size(p, 2) != 3
            error("Expected points matrix with 3 columns, got size $(size(p))")
        end
        
        distances = zeros(size(p, 1))
        p1 = [x1, y1, z1]
        
        for i in 1:size(p, 1)
            pt = p[i, :]
            v = pt - p1
            proj = dot(v, axis)
            proj = clamp(proj, 0, len)
            closest = p1 + proj * axis
            radial_dist = norm(pt - closest)
            distances[i] = radial_dist - r
        end
        
        return distances
    end
    
    return (surface, sdf)
end

end  # module