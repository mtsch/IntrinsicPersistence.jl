using IntrinsicPersistence
using Test

using StaticArrays

# sanitycheck
@testset "Sanity checks" begin
    @testset "regular circle" begin
        n = 12
        pts = [SVector(cos(t), sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]]
        landmarks = 1:2:n

        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test length(res.cycles) == 1
        α = first(res.cycles).points
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
        @test length(res.cycles) == 1
        α = first(res.cycles).points
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
        @test length(res.cycles) == 0
    end

    @testset "cylinder" begin
        n = 12
        pts = [SVector(cos(t), sin(t), h)
               for t in range(0, stop = 2π, length = n+1)[1:end-1]
               for h in 0:0.5:10]
        landmarks = 1:2:length(pts)
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test length(res.cycles) == 1
    end

    @testset "???" begin
        n = 12
        pts = [SVector(h/10*cos(t), h/10*sin(t), h)
               for t in range(0, stop = 2π, length = n+1)[1:end-1]
               for h in 0:0.5:10]
        landmarks = 1:2:length(pts)
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test_broken length(res.cycles) == 0
    end
end
