using CollocationPoints
using Test
using LinearAlgebra

@testset "CollocationPoints.jl" begin
    
    @testset "draw_circ" begin
        xc, yc, r, n = 0.5, -0.5, 2.0, 100
        xy, sdf = draw_circ(xc, yc, r, n)
        
        @test size(xy, 1) == n
        @test size(xy, 2) == 2
        
        # Test if points are on the circle
        distances = sqrt.((xy[:,1] .- xc).^2 .+ (xy[:,2] .- yc).^2)
        @test all(d -> isapprox(d, r, atol=1e-9), distances)
        
        # Test signed distance function
        @test isapprox(sdf([xc yc]), -r)
        @test isapprox(sdf([xc+r yc]), 0.0)
    end

    @testset "draw_rect" begin
        x1, y1, x2, y2, n = -1.0, -1.0, 1.0, 1.0, 10
        bdy, sdf = draw_rect(x1, y1, x2, y2, n)
        
        @test size(bdy, 1) == 4 * (n-1)
        @test size(bdy, 2) == 2
        
        # Test signed distance function
        @test sdf([0.0 0.0]) < 0 # inside
        @test sdf([2.0 2.0]) > 0 # outside
        @test isapprox(sdf([1.0 0.0]), 0.0) # on boundary
    end

    @testset "bsmooth" begin
        # Test with a simple square
        square = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.0 0.0]
        smoothed = bsmooth(square, 0.1)
        @test size(smoothed, 1) > size(square, 1)
        
        # Test it doesn't error on small inputs
        @test size(bsmooth([0.0 0.0], 0.1), 1) == 1
        @test size(bsmooth(zeros(0,2), 0.1), 1) == 0
    end

end