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
        n = 24
        pts = vcat(
            [SVector(1.0cos(t), 1.0sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]],
            [SVector(0.9cos(t), 0.9sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]],
            [SVector(0.8cos(t), 0.8sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]],
            [SVector(0.7cos(t), 0.7sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]],
            [SVector(0.6cos(t), 0.6sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]],
            [SVector(0.5cos(t), 0.5sin(t)) for t in range(0, stop = 2π, length = n+1)[1:end-1]],
        )
        landmarks = 1:2:n
        gc = GeodesicComplex(pts, 0.3, landmarks = landmarks)
        res = persistence(gc)
        @test length(res) == 1
        @test diameter(first(res)) ≈ 1
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

        res = persistence(gc, shrinked = false)
        @test !isempty(res)
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

    @testset "random full graph" begin
        n = 100
        pts = [SVector(rand(), rand(), rand()) for _ in 1:n]
        landmarks = 1:n
        gc = GeodesicComplex(pts, 2, landmarks = landmarks)
        res = persistence(gc)
        @test isempty(res)
    end

    @testset "regular cone" begin
        n = 12
        pts = [SVector(h/10*cos(t), h/10*sin(t), h)
               for t in range(0, stop = 2π, length = n+1)[1:end-1]
               for h in 0:0.5:10]
        landmarks = 1:2:length(pts)
        gc = GeodesicComplex(pts, 0.6, landmarks = landmarks)
        res = persistence(gc)
        @test isempty(res)
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
