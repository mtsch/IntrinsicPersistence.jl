"""
    getindexmatrix(dists)

Create index matrix `M`, such that `M[i, j] = k` if the value `dists[i, j]` is the `k`-th
largest. Assumes `dists` is symmetric.
"""
function getindexmatrix(dists)
    M = fill(0, size(dists))
    for (i, (u, v)) in enumerate(sortededges(dists))
        M[u, v] = i
        M[v, u] = i
    end
    M
end

"""
Internal state for intrinsic persistence.
"""
struct PersistenceState{T}
    # Precomputed properties.
    dists       ::Matrix{T}
    parents     ::Matrix{Int}
    indexmatrix ::Matrix{Int}
    triangles   ::EquilateralIterator{T}
    radius      ::T
    # Internal state
    reduced     ::Vector{BitSet}
    # Buffers
    σ           ::BitSet
    cycle       ::Vector{Int}
    # Results
    results     ::Vector{Tuple{Vector{Int}, NTuple{3, Int}, T}}
end

function PersistenceState(gc)
    # todo: is this correct? should be maximum(weights)?
    r = radius(gc)
    # todo: split connected components
    dists, parents = landmark_shortest_paths(gc)
    indexmatrix = getindexmatrix(dists)
    triangles = equilaterals(dists, 2r)

    reduced = fill(BitSet(), binomial(n_landmarks(gc), 2))

    PersistenceState(dists, parents, indexmatrix, triangles, r, reduced,
                     BitSet(), Int[], Tuple{Vector{Int}, NTuple{3, Int}, eltype(dists)}[])
end

Base.show(io::IO, st::PersistenceState) = print(io, "PersistenceState")

const NOP_COUNT = Ref(0)
const YEP_COUNT = Ref(0)
const ALL_COUNT = Ref(0)
const RED_COUNT = Ref(0)

"""
    isgeodesic(st, Δ)

Check if `st.cycle` spanned by triangle `Δ` is a geodesic circle.
"""
function isgeodesic(st, Δ)
    length(st.cycle) == 3 && return true
    d = st.dists
    dik = djk = zero(eltype(d))
    ALL_COUNT[] += 1
    # add offset to get it right in the first go.
    j, k, i = Δ
    for x in st.cycle
        if x ∈ Δ
            # make sure x is always between i and j.
            i, j, k = j, k, i
            dik = d[i, k]
            djk = d[j, k]
            continue
        end
        # add st.radius/2 to a void numerical errors.
        # st.radius is the shortest edge in graph.
        if d[x, k] + st.radius/2 < min(d[x, i] + dik,
                                       d[x, j] + djk)
            NOP_COUNT[] += 1
            return false
        end
    end
    YEP_COUNT[] += 1
    true
end

"""
    low(σ::BitSet)

Get the largest index in `σ` / the lowest entry in the column.
"""
low(σ::BitSet) = length(σ) > 0 ? last(σ) : 0

"""
    reduce!(st::PersistenceState)

Use the standard reduction algorithm to reduce `st.σ`.
"""
function reduce!(st::PersistenceState)
    l = low(st.σ)
    while l ≠ 0 && !isempty(st.reduced[l])
        symdiff!(st.σ, st.reduced[l])
        l = low(st.σ)
    end
    l
end

"""
    appendpath!(st::PersistenceState, u, v)

Append the shortest path from `u` to `v` into `st.cycle` and `st.σ`.
"""
function appendpath!(st::PersistenceState, u, v)
    while u ≠ v
        u′ = st.parents[v, u]
        push!(st.σ, st.indexmatrix[u, u′])
        push!(st.cycle, u)
        u = u′
    end
end

"""
    movecycle!(st::PersistenceState, Δ)

Move the cycle spanned by triangle `Δ` into `st.cycle` and `st.σ`.
"""
function movecycle!(st::PersistenceState, (i, j, k))
    empty!(st.σ)
    empty!(st.cycle)
    appendpath!(st, i, j)
    appendpath!(st, j, k)
    appendpath!(st, k, i)
end

"""
TODO
"""
function processtriangle!(st::PersistenceState, triangle, diameter)
    movecycle!(st, triangle)
    # Skip cycles that visit a node more than once and nongeodesic cycles.
    length(st.σ) == length(st.cycle) || return
    isgeodesic(st, triangle) || return

    l = reduce!(st)
    if l ≠ 0
        st.reduced[l] = copy(st.σ)
        push!(st.results, (copy(st.cycle), triangle, diameter))
    else
        RED_COUNT[] += 1
    end
end

function persistence(gc; showprogress = false, dense = false)
    NOP_COUNT[] = 0
    YEP_COUNT[] = 0
    ALL_COUNT[] = 0
    RED_COUNT[] = 0
    showprogress && println("Calculating intrinsic persistence...")
    st = PersistenceState(gc)

    for (Δ, diam) in st.triangles
        processtriangle!(st, Δ, diam)
    end
    showprogress && println("Postprocessing...")
    cycles = map(st.results) do (α, _, _)
        landmarks(gc, α)
    end
    postprocess(gc, cycles, dense = dense)
end

# ======================================================================================== #

"""
find the index and distance of point in `gc` within `gc.radius` of point(s) `p` that
minimizes the maximum distance to points in `pts`.
"""
function centerpoint(gc, p, pts)
    inball = nearby_points(gc, p)
    dists = pairwise_ambient_distance(gc, inball, pts)
    μ, i = findmin(vec(maximum(dists, dims = 2)))
    μ, inball[i]
end

shrink(gc, α) where {T} = shrink!(gc, copy(α))

function shrink!(gc, α)
    changed = true
    while changed
        changed = false
        unique!(α)
        n = length(α)
        for (i, v) in enumerate(α)
            j = mod1(i - 1, n)
            k = mod1(i + 1, n)
            _, new = centerpoint(gc, v, [α[j], α[k]])
            if new ≠ v
                changed = true
                α[i] = new
            end
        end
    end
    α
end

# Kaj bo s tem?
# TODO: densify + shrink
function densify(gc::GeodesicComplex, α)
    α = copy(α)
    β = Int[]
    changed = true
    while changed
        changed = false
        for i in eachindex(α)
            u = α[i]
            v = α[mod1(i + 1, length(α))]
            push!(β, u)
            inballs = nearby_points(gc, [u, v])
            candidates = union(inballs...)
            dists = pairwise_ambient_distance(gc, candidates, [u, v])
            new = candidates[argmin(vec(maximum(dists, dims = 2)))]
            if new ≠ u && new ≠ v
                changed = true
                push!(β, new)
            end
        end
        move!(α, β)
        shrink!(gc, α)
    end
    α
end

#=
function criticalpoints(gc::GeodesicComplex, α)
    α = copy(α)
    ε = zero(T)
    changed = true
    while changed
        changed = false
        d = ε = typemax(T)
        for (i, u) in enumerate(α)
            inball = nearby_points(gc, u)
            # Spaghetti mapreduce avoids allocating distance matrices.
            ε, j = mapreduce(min, enumerate(inball)) do (i, v)
                (maximum(evaluate(gc.metric, points(gc, v), landmarks(gc, w))
                         for w in cycle), i)
            end
            d = min(ε, d)
            changed |= α[i] ≠ N[j]
            α[i] = N[j]
        end
        unique!(α)
    end
    α, ε
end
=#

function move!(a, b)
    resize!(a, length(b))
    copyto!(a, b)
    empty!(b)
    a
end

function postprocess(gc, cycles; dense = false)
    res = map(cycles) do α
        if dense
            cycle = densify(gc, α)
        else
            cycle = shrink(gc, α)
        end
        Cycle(gc, cycle, perimeter(gc, cycle), diameter(gc, cycle))
    end
    filter!(res) do cycle
        length(cycle) > 2 && diameter(cycle) > radius(gc)
    end
end

#TODO new interface
perimeter(gc, α) =
    sum(eachindex(α)) do i
        evaluate(gc.metric, points(gc, α[i]), points(gc, α[mod1(i + 1, length(α))]))
    end

LightGraphs.diameter(gc, α) = maximum(pairwise_ambient_distance(gc, α, α))

struct Cycle{T, G<:GeodesicComplex{T}}
    complex   ::G
    indices   ::Vector{Int}
    perimeter ::T
    diameter  ::T
end

Base.show(io::IO, α::Cycle) =
    print(io, "$(length(α))-element Cycle with p = $(perimeter(α)), d = $(diameter(α))")

Base.length(α::Cycle) = length(α.indices)

points(α::Cycle) = points(α.complex, α.indices)

perimeter(α::Cycle) = α.perimeter

LightGraphs.diameter(α::Cycle) = α.diameter #?

#TODO densify(cycle, niter=Inf)
# TODO: critical points
