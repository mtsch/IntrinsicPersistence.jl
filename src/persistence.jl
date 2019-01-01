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

Get the largest index in `σ`.
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
    movepath!(st::PersistenceState, u, v)

Move shortest path from `u` to `v` into `st.cycle` and `st.σ`.
"""
function movepath!(st::PersistenceState, u, v)
    while u ≠ v
        u′ = st.parents[v, u]
        push!(st.σ, st.indexmatrix[u, u′])
        push!(st.cycle, u)
        u = u′
    end
end

"""
    movecycle!(st::PersistenceState, Δ)

Move cycle spanned by `Δ` into `st.cycle` and `st.σ`.
"""
function movecycle!(st::PersistenceState, (i, j, k))
    empty!(st.σ)
    empty!(st.cycle)
    movepath!(st, i, j)
    movepath!(st, j, k)
    movepath!(st, k, i)
end

"""
TODO
"""
function processtriangle!(st::PersistenceState, triangle, diameter, showprogress)
    movecycle!(st, triangle)
    # Skip cycles that visit a node more than once.
    length(st.σ) == length(st.cycle) || return
    # Geodesic check
    isgeodesic(st, triangle) || return

    l = reduce!(st)
    if l ≠ 0
        st.reduced[l] = copy(st.σ)
        push!(st.results, (copy(st.cycle), triangle, diameter))
    else
        #println(st.cycle)
        RED_COUNT[] += 1
    end
end

function persistence(gc, showprogress = false)
    NOP_COUNT[] = 0
    YEP_COUNT[] = 0
    ALL_COUNT[] = 0
    RED_COUNT[] = 0
    showprogress && println("Calculating intrinsic persistence...")
    st = PersistenceState(gc)

    for (Δ, diam) in st.triangles
        processtriangle!(st, Δ, diam, showprogress)
    end
    showprogress && println("Postprocessing...")
    cycles = map(st.results) do (α, Δ, d)
        landmarks(gc, α), landmarks.(Ref(gc), Δ), d
    end
    postprocess(gc, cycles)
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

function shrink(gc::GeodesicComplex{T}, α) where {T}
    α = copy(α)
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
    α, maximum(pairwise_ambient_distance(gc, α, α))
end

# Kaj bo s tem?
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
    end
    α
end

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

function move!(a, b)
    resize!(a, length(b))
    copyto!(a, b)
    empty!(b)
    a
end

function postprocess(gc, cycles)
    res = map(cycles) do (α, Δ, d)
        cycle, diam = shrink(gc, α)
        #dense = densify(gc, cycle)
        dense = copy(α)
        Cycle(cycle, dense, Δ, perimeter(gc, cycle), diam, d, 0)
    end
    filter!(res) do cycle
        length(cycle.points) > 2 || cycle.diameter > gc.radius
    end
    IntrinsicPersistenceResults(gc, res)
end

function perimeter(gc::GeodesicComplex{T}, α) where {T}
    res = zero(T)
    for i in eachindex(α)
        res += evaluate(gc.metric, points(gc, α[i]), points(gc, α[mod1(i + 1, length(α))]))
    end
    res
end

struct Cycle{T}
    points     ::Vector{Int}
    dense      ::Vector{Int}
    triangle   ::NTuple{3, Int}
    perimeter  ::T
    diameter   ::T
    death      ::T
    deathpoint ::Int
end

struct IntrinsicPersistenceResults{T, G<:GeodesicComplex{T}}
    complex    ::G
    cycles ::Vector{Cycle{T}}
end

# TODO: critical points
