using IntrinsicPersistence
using Test

using StaticArrays
using LinearAlgebra

# sanitycheck
@testset "Sanity checks" begin
    @testset "regular circle" begin
        n = 12
        pts = [SVector(cos(t), sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]]
        landmarks = 1:2:n

        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test length(res) == 1
        α = first(res).indices
        @test sort(α) == landmarks
    end

    @testset "circle with contraction" begin
        n = 12
        pts = vcat([SVector(0.5cos(t), 0.5sin(t))
                    for t in range(0, stop = 2π, length = n÷2+1)[1:end-1]],
                   [SVector(cos(t), sin(t))
                    for t in range(0, stop = 2π, length = n+1)[1:end-1]])
        landmarks = n÷2+1:2:n+n÷2
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test length(res) == 1
        α = first(res).indices
        @test sort(α) == 1:n÷2
    end

    @testset "landmarks form a circle, but the underlying space is contractible" begin
        n = 12
        pts = vcat([SVector(0.5cos(t), 0.5sin(t))
                    for t in range(0, stop = 2π, length = n÷2+1)[1:end-1]],
                   [SVector(cos(t), sin(t))
                    for t in range(0, stop = 2π, length = n+1)[1:end-1]],
                   [SVector(0.0, 0.0)])
        landmarks = n÷2+1:2:n+n÷2
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test isempty(res)
    end

    @testset "cylinder" begin
        n = 12
        pts = [SVector(cos(t), sin(t), h)
               for t in range(0, stop = 2π, length = n+1)[1:end-1]
               for h in 0:0.5:10]
        landmarks = 1:2:length(pts)
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test length(res) == 1
    end

    #=
    @testset "cone, landmarks in a circle" begin
        n = 12
        pts = [SVector((1 - h/10)cos(t), (1 - h/10)sin(t), h)
               for h in 0:0.5:10
               for t in range(0, stop = 2π, length = n+1)[1:end-1]]
        landmarks = 1:2:n
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test_broken length(res.cycles) == 0
    end

    @testset "truncated cone, landmarks in a circle" begin
        n = 12
        pts = [SVector((2 - h/10)*cos(t), (2 - h/10)*sin(t), h)
               for h in 0:0.5:10
               for t in range(0, stop = 2π, length = n+1)[1:end-1]]
        landmarks = 1:2:n
        gc = GeodesicComplex(pts, 1.2, landmarks = landmarks)
        res = persistence(gc)
        @test length(res.cycles) == 1
        # diameter = 2
    end
    =#

    @testset "random full graph" begin
        n = 100
        pts = [SVector(rand(), rand(), rand()) for _ in 1:n]
        landmarks = 1:n
        gc = GeodesicComplex(pts, 2, landmarks = landmarks)
        res = persistence(gc)
        @test isempty(res)
    end

    @testset "???" begin
        n = 12
        pts = [SVector(h/10*cos(t), h/10*sin(t), h)
               for t in range(0, stop = 2π, length = n+1)[1:end-1]
               for h in 0:0.5:10]
        landmarks = 1:2:length(pts)
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test_broken isempty(res)
    end

    @testset "random cone" begin
        n = 100_000
        pts = map(rand() for _ in 1:n) do h
            t = rand()
            SVector(h*sinpi(2t), h*cospi(2t), h)
        end
        gc = GeodesicComplex(pts, 0.2)
        res = persistence(gc)
        @test isempty(res)
    end

    @testset "random truncated cone" begin
        n = 100_000
        pts = map(rand() for _ in 1:n) do h
            t = rand()
            SVector((h+1)*sinpi(2t), (h+1)*cospi(2t), h)
        end
        gc = GeodesicComplex(pts, 0.2)
        res = persistence(gc)
        @test length(res) == 1
        @test diameter(first(res)) ≈ 2 atol = 0.2
        @test perimeter(first(res)) ≈ 2π atol = 0.2
    end
end
