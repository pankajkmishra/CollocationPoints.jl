module CollocationPoints

using LinearAlgebra
using Random
using DelimitedFiles
using Interpolations

export CollocationPoints2D, draw_circ, draw_rect, make_domain, bsmooth, min_pdist2

function node_placing(box, ninit, dotmax, ctps, radius_fn)
    dotnr = 0
    Random.seed!(0)
    pdp_x = range(box[1], stop=box[2], length=ninit)
    pdp_y = box[3] .+ 1e-4 .* rand(ninit)
    pdp = hcat(pdp_x, pdp_y)
    p = zeros(dotmax, 2)
    if isempty(pdp)
        return zeros(0, 2)
    end
    ym, i = findmin(pdp[:, 2])
    while ym <= box[4] && dotnr < dotmax
        dotnr += 1
        p[dotnr, :] = pdp[i, :]
        current_point = reshape(p[dotnr, :], 1, 2)
        r = radius_fn(current_point, ctps)[1]
        dist2 = sum((pdp .- p[dotnr, :]').^2, dims=2)[:]
        ileft_candidates = findall(dist2[1:i] .> r^2)
        ileft = isempty(ileft_candidates) ? 0 : last(ileft_candidates)
        ang_left = (ileft == 0) ? pi : atan(pdp[ileft, 2] - pdp[i, 2], pdp[ileft, 1] - pdp[i, 1])
        iright_candidates = findall(@view(dist2[i:end]) .> r^2)
        iright_relative = isempty(iright_candidates) ? 0 : first(iright_candidates)
        actual_iright_index = i + iright_relative - 1
        ang_right = (iright_relative == 0) ? 0.0 : atan(pdp[actual_iright_index, 2] - pdp[i, 2], pdp[actual_iright_index, 1] - pdp[i, 1])
        ang = ang_left .- [0.1, 0.3, 0.5, 0.7, 0.9] .* (ang_left - ang_right)
        pdp_new = hcat(pdp[i, 1] .+ r .* cos.(ang), pdp[i, 2] .+ r .* sin.(ang))
        ind_to_keep = (pdp_new[:, 1] .>= box[1]) .& (pdp_new[:, 1] .<= box[2])
        pdp_new = pdp_new[ind_to_keep, :]
        pdp_top = pdp[1:ileft, :]
        pdp = (iright_relative == 0) ? vcat(pdp_top, pdp_new) : vcat(pdp_top, pdp_new, pdp[actual_iright_index:end, :])
        if isempty(pdp)
            break
        end
        ym, i = findmin(pdp[:, 2])
    end
    return p[1:dotnr, :]
end

function CollocationPoints2D(fd, bbox, ctps, ptol, radius_fn)
    ninit = 50000
    dotmax = 500000
    p = node_placing(bbox, ninit, dotmax, ctps, radius_fn)
    if isempty(p) return zeros(0,2) end
    p = p[fd(p) .< eps(), :]
    if !isempty(ctps) && !isempty(p)
        p_tuples = [tuple(row...) for row in eachrow(p)]
        ctps_tuples = [tuple(row...) for row in eachrow(ctps)]
        p_tuples = setdiff(p_tuples, ctps_tuples)
        p = isempty(p_tuples) ? zeros(0, 2) : vcat([collect(t)' for t in p_tuples]...)
    end
    pfix = isempty(ctps) ? zeros(0,2) : unique(ctps, dims=1)
    p = vcat(pfix, p)
    if isempty(p) return zeros(0,2) end
    ib = findall(abs.(fd(p)) .< ptol)
    p = p[setdiff(1:size(p,1), ib), :]
    return p
end

function dsegment(p, pv)
    np = size(p, 1)
    nvs = size(pv, 1)
    ds = zeros(np, nvs - 1)
    for iv in 1:(nvs-1)
        v1 = pv[iv, :]
        v2 = pv[iv+1, :]
        v = v2 - v1
        c2 = dot(v, v)
        if c2 < 1e-12
            for ip in 1:np
                ds[ip, iv] = norm(p[ip, :] - v1)
            end
            continue
        end
        for ip in 1:np
            w = p[ip, :] - v1
            c1 = dot(v, w)
            if c1 <= 0
                ds[ip, iv] = norm(p[ip, :] - v1)
            elseif c1 >= c2
                ds[ip, iv] = norm(p[ip, :] - v2)
            else
                pb = v1 + (c1 / c2) * v
                ds[ip, iv] = norm(p[ip, :] - pb)
            end
        end
    end
    return ds
end

function inpolygon(p, poly)
    n = size(poly, 1)
    inside = falses(size(p, 1))
    for i in 1:size(p, 1)
        point = p[i,:]
        x, y = point[1], point[2]
        is_inside = false
        p1x, p1y = poly[n, :]
        for j in 1:n
            p2x, p2y = poly[j, :]
            if y > min(p1y, p2y) && y <= max(p1y, p2y) && x <= max(p1x, p2x)
                if p1y != p2y
                    xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x || x <= xinters
                        is_inside = !is_inside
                    end
                elseif x <= p1x
                     is_inside = !is_inside
                end
            end
            p1x, p1y = p2x, p2y
        end
        inside[i] = is_inside
    end
    return inside
end

function getsdf(xy, bdy)
    d = dsegment(xy, bdy)
    s = (-1) .^ inpolygon(xy, bdy)
    return s .* minimum(d, dims=2)[:]
end

function draw_circ(xc, yc, r, n)
    n_int = round(Int, n)
    th = range(0, 2*pi, length=n_int+1)[1:n_int]
    xy = hcat(r .* cos.(th) .+ xc, r .* sin.(th) .+ yc)
    sdf(p) = @. sqrt((p[:,1] - xc)^2 + (p[:,2] - yc)^2) - r
    return (xy, sdf)
end

function draw_rect(x1, y1, x2, y2, n)
    n_int = round(Int, n)
    b1 = hcat(range(x1, x2, length=n_int), fill(y1, n_int))
    b2 = hcat(fill(x2, n_int), range(y1, y2, length=n_int))
    b3 = hcat(range(x2, x1, length=n_int), fill(y2, n_int))
    b4 = hcat(fill(x1, n_int), range(y2, y1, length=n_int))
    B = vcat(b1, b2[2:end,:], b3[2:end,:], b4[2:end-1,:])
    sdf(p) = @. -min(min(min(-y1 + p[:,2], y2 - p[:,2]), -x1 + p[:,1]), x2 - p[:,1])
    return (B, sdf)
end

function make_domain(fileName)
    xy = readdlm(fileName)
    sdf(p) = getsdf(p, xy)
    return (xy, sdf)
end

function min_pdist2(points, ctps)
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

# Simple implementation of bsmooth that doesn't rely on complex interpolation methods
function bsmooth(b, h0)
    if size(b, 1) < 2
        return b
    end
    
    # Remove duplicate points
    unique_indices = [1]
    for i in 2:size(b, 1)
        if !isapprox(b[i,:], b[unique_indices[end],:])
            push!(unique_indices, i)
        end
    end
    b_clean = b[unique_indices, :]

    if size(b_clean, 1) < 2
        return b_clean
    end

    # Check if boundary is closed
    is_closed = isapprox(b_clean[1,:], b_clean[end,:])

    if is_closed
        # For closed boundaries, use only points up to the last one (which is a duplicate of the first)
        b_interp = b_clean[1:end-1, :]
    else
        b_interp = b_clean
    end
    
    if size(b_interp, 1) < 2
        return b_clean
    end

    # Calculate cumulative distances along the boundary
    t = vcat(0, cumsum(norm.(eachrow(diff(b_interp, dims=1)))))
    
    # Calculate total boundary length
    if is_closed
        total_length = last(t) + norm(b_interp[end,:] - b_interp[1,:])
    else
        total_length = last(t)
    end
    
    if total_length < 1e-9 
        return b_interp
    end

    # Calculate number of points for the smoothed boundary
    num_points = max(100, round(Int, total_length / (h0/2)))
    
    # Generate new points along the boundary using linear interpolation
    if is_closed
        # For closed curve, wrap around
        t_new = range(0, stop=total_length, length=num_points+1)[1:end-1]
        b_new = zeros(num_points, 2)
        
        for i in 1:num_points
            t_i = t_new[i]
            while t_i > last(t)
                t_i -= last(t)
            end
            
            # Find segment containing t_i
            seg_idx = findlast(t .<= t_i)
            if seg_idx === nothing || seg_idx >= length(t)
                seg_idx = length(b_interp) - 1
                next_idx = 1
            else
                next_idx = seg_idx + 1
                if next_idx > length(b_interp)
                    next_idx = 1
                end
            end
            
            # Linear interpolation
            if seg_idx < length(t)
                alpha = (t_i - t[seg_idx]) / (t[seg_idx+1] - t[seg_idx])
                b_new[i,:] = (1 - alpha) * b_interp[seg_idx,:] + alpha * b_interp[next_idx,:]
            else
                # Handle wrapping around to the beginning
                alpha = (t_i - t[end]) / (total_length - t[end])
                b_new[i,:] = (1 - alpha) * b_interp[end,:] + alpha * b_interp[1,:]
            end
        end
    else
        # For open curve
        t_new = range(0, stop=total_length, length=num_points)
        b_new = zeros(num_points, 2)
        
        for i in 1:num_points
            t_i = t_new[i]
            if t_i >= last(t)
                b_new[i,:] = b_interp[end,:]
                continue
            end
            
            # Find segment containing t_i
            seg_idx = findlast(t .<= t_i)
            if seg_idx === nothing
                seg_idx = 1
            end
            next_idx = min(seg_idx + 1, length(b_interp))
            
            # Linear interpolation
            if seg_idx < length(t)
                alpha = (t_i - t[seg_idx]) / (t[seg_idx+1] - t[seg_idx])
                b_new[i,:] = (1 - alpha) * b_interp[seg_idx,:] + alpha * b_interp[next_idx,:]
            else
                b_new[i,:] = b_interp[end,:]
            end
        end
    end
    
    return b_new
end

end