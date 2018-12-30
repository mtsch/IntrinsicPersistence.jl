module GeodesicComplexes
export GeodesicComplex, npoints, nlandmarks, points, landmarks, nonlandmarks, landmark_idxs

using Distances
using LightGraphs
using NearestNeighbors
using RecipesBase
using SimpleWeightedGraphs
using StaticArrays

using Random
using SparseArrays

struct GeodesicComplex{T, P<:AbstractVector, M<:Metric, K<:NNTree{P, M}} <:
                       AbstractSimpleWeightedGraph{Int, T}
    points        ::Vector{P}
    landmark_idxs ::Vector{Int}
    graph         ::SimpleWeightedGraph{Int, T} # <- rab se samo dists in paths
    tree          ::K
    cover         ::Vector{Set{Int}} # <- se ne rab
    radius        ::T
    metric        ::M
end

function GeodesicComplex(pts::AbstractVector{<:AbstractVector}, r; landmarks = nothing,
                         metric = Euclidean(), tree = KDTree)
    kdt = tree(pts, metric)
    landmark_idxs, cover = getcover(landmarks, pts, r, kdt)
    I = Int[]
    J = Int[]
    T = result_type(metric, pts[1], pts[1])
    V = T[]
    n = length(landmark_idxs)

    for i in 1:n, j in 1:(i - 1)
        if !isempty(intersect(cover[i], cover[j]))
            d = evaluate(metric, pts[landmark_idxs[i]], pts[landmark_idxs[j]])
            append!(I, (i, j))
            append!(J, (j, i))
            append!(V, (d, d))
        end
    end
    graph = SimpleWeightedGraph(sparse(I, J, V, n, n))

    GeodesicComplex(pts, landmark_idxs, graph, kdt, cover, T(r), metric)
end

function getcover(::Nothing, pts, r, tree)
    covered = falses(length(pts))
    idxs = shuffle(eachindex(pts))
    landmarks = Int[]
    cover = Set{Int}[]

    for i in idxs
        covered[i] && continue
        covered[i] = true
        inball = inrange(tree, pts[i], r)
        covered[inball] .= true
        push!(landmarks, i)
        push!(cover, Set(inball))
    end
    ord = sortperm(landmarks)
    landmarks[ord], cover[ord]
end

function getcover(landmarks, pts, r, tree)
    cover = Set{Int}[]

    for l in landmarks
        inball = inrange(tree, pts[l], r)
        push!(cover, Set(inball))
    end
    ord = sortperm(landmarks)
    landmarks[ord], cover[ord]
end

# Getters
# TODO = Int?
Base.eltype(::GeodesicComplex{<:Any, P}) where {P} = P

npoints(gc::GeodesicComplex) =
    length(gc.points)

nlandmarks(gc::GeodesicComplex) =
    length(gc.landmark_idxs)

points(gc::GeodesicComplex) =
    gc.points
points(gc::GeodesicComplex, idxs) =
    gc.points[idxs]

landmarks(gc::GeodesicComplex, idxs=1:nlandmarks(gc)) =
    gc.points[gc.landmark_idxs[idxs]]

landmark_idxs(gc::GeodesicComplex) =
    gc.landmark_idxs
landmark_idxs(gc::GeodesicComplex, idxs) =
    gc.landmark_idxs[idxs]

nonlandmarks(gc::GeodesicComplex, idxs=1:npoints(gc)-nlandmarks(gc)) =
    gc.points[setdiff(1:npoints(gc), gc.landmark_idxs)[idxs]]

points_mat(gc::GeodesicComplex{<:Any, P}, idxs) where {P} =
    resape(reinterpret(T, points(gc, idxs)), (length(P), length(idxs)))

# LightGraphs interface
for f in [:edges, :ne, :nv, :vertices, :is_directed, :weights]
    @eval LightGraphs.$f(gc::GeodesicComplex) = $f(gc.graph)
end
LightGraphs.is_directed(::Type{GeodesicComplex}) = false
LightGraphs.is_directed(::Type{GeodesicComplex{T,P,M,K}}) where {T,P,M,K} = false
LightGraphs.has_edge(gc::GeodesicComplex, u::Integer, v::Integer) = has_edge(gc.graph, u, v)
LightGraphs.has_vertex(gc::GeodesicComplex, v::Integer) = has_edge(gc.graph, v)
LightGraphs.inneighbors(gc::GeodesicComplex, v::Integer) = inneighbors(gc.graph, v)
LightGraphs.outneighbors(gc::GeodesicComplex, v::Integer) = outneighbors(gc.graph, v)

function Base.show(io::IO, gc::GeodesicComplex{T, P}) where {T, P}
    print(io, "GeodesicComplex{$T, $P} with $(npoints(gc)) points, " *
          "$(nlandmarks(gc)) landmarks and $(ne(gc)) edges")
end

# NearestNeighbors interface
const NN = NearestNeighbors
NN.inrange(gc::GeodesicComplex, i, r = gc.radius, sortres = false) =
    inrange(gc.tree, i, r, sortres)
NN.knn(gc::GeodesicComplex, pts, k, sortres = false, skip = NN.always_false) =
    knn(gc.tree, pts, k, sortres, skip)

# DistancesInterface
Distances.result_type(::GeodesicComplex{T}, ::Int, ::Int) where {T} = T
Distances.evaluate(gc::GeodesicComplex, i, j) =
    evaluate(gc.metric, points(gc, i), points(gc, j))
Distances.colwise!(res, gc::GeodesicComplex, is, js) =
    colwise!(res, gc.metric, is, js)
Distances.colwise(gc::GeodesicComplex, is, js) =
    colwise(gc.metric, is, js)
Distances.pairwise!(res, gc::GeodesicComplex, is, js) =
    pairwise!(res, gc.metric, is, js)
Distances.pairwise(gc::GeodesicComplex, is, js) =
    pairwise(gc.metric, is, js)

# Plots recipe
function getxyz(pts)
    xs = get.(pts, 1, 0.0)
    ys = get.(pts, 2, 0.0)
    zs = get.(pts, 3, 0.0)
    if all(iszero, zs)
        xs, ys
    else
        xs, ys, zs
    end
end

@recipe function plot(gc::GeodesicComplex; only_landmarks = true, graph = true)
    # edges
    if graph
        @series begin
            label := "edges"
            edgepoints = Vector{Float64}[]
            for e in edges(gc)
                s = landmarks(gc, src(e))
                d = landmarks(gc, dst(e))
                append!(edgepoints, (s, d, fill(NaN, length(s))))
            end
            getxyz(edgepoints)
        end
    end
    # landmarks
    @series begin
        markersize --> 1.0
        seriestype --> :scatter
        label := "landmarks"
        getxyz(landmarks(gc))
    end
    # others
    if !only_landmarks
        @series begin
            markersize --> 0.5
            seriestype --> :scatter
            label := "others"
            getxyz(points(gc, setdiff(1:npoints(gc), landmark_idxs(gc))))
        end
    end
end

end
