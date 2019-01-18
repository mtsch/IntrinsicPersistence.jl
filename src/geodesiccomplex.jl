using Distances
using LightGraphs
using NearestNeighbors
using SimpleWeightedGraphs
using StaticArrays

using Random
using SparseArrays

struct GeodesicComplex{T, P<:AbstractVector, M<:Metric,
                       I<:AbstractVector{Int}, K<:NNTree{P, M}}
    points    ::Vector{P}
    radius    ::T
    metric    ::M
    landmarks ::I
    graph     ::SimpleWeightedGraph{Int, T} # <- rab se samo dists in paths
    tree      ::K
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

    GeodesicComplex(pts, T(r), metric, landmark_idxs, graph, kdt)
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
    landmarks = sort(landmarks)

    for l in landmarks
        inball = inrange(tree, pts[l], r)
        push!(cover, Set(inball))
    end
    landmarks, cover
end

function Base.show(io::IO, gc::GeodesicComplex{T, P}) where {T, P}
    print(io, "GeodesicComplex{$T, $P} with $(n_points(gc)) points, " *
          "$(n_landmarks(gc)) landmarks and $(ne(landmark_graph(gc))) edges")
end

# new interface:
landmarks(gc::GeodesicComplex) = gc.landmarks
landmarks(gc::GeodesicComplex, idxs) = gc.landmarks[idxs]

n_landmarks(gc::GeodesicComplex) = length(gc.landmarks)

points(gc::GeodesicComplex) = gc.points
points(gc::GeodesicComplex, idxs) = gc.points[idxs]

n_points(gc::GeodesicComplex) = length(gc.points)

radius(gc::GeodesicComplex) = gc.radius

landmark_graph(gc::GeodesicComplex) = gc.graph

function landmark_shortest_paths(gc)
    fw = floyd_warshall_shortest_paths(landmark_graph(gc))
    fw.dists, fw.parents
end

pairwise_ambient_distance(gc::GeodesicComplex{T, D}, is, js) where {T, D} =
    pairwise(gc.metric,
             reshape(reinterpret(T, points(gc, is)), (length(D), length(is))),
             reshape(reinterpret(T, points(gc, js)), (length(D), length(js))))

nearby_points(gc::GeodesicComplex, i, r = radius(gc)) = inrange(gc.tree, points(gc, i), r)

distance_result_type(gc::GeodesicComplex{T}) where T = T
