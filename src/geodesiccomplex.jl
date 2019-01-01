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


#=
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
    gc.points[setdiff(1:npoints(gc), gc.landmark)[idxs]]

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

# NearestNeighbors interface
const NN = NearestNeighbors
NN.inrange(gc::GeodesicComplex, i, r = gc.radius, sortres = false) =
    inrange(gc.tree, i, r, sortres)
NN.knn(gc::GeodesicComplex, pts, k, sortres = false, skip = NN.always_false) =
    knn(gc.tree, pts, k, sortres, skip)

# Distances interface
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

=#
