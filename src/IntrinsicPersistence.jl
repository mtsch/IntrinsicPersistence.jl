module IntrinsicPersistence

using Distances
using LightGraphs
using NearestNeighbors
using RecipesBase
using SimpleWeightedGraphs
using StaticArrays

using Random
using SparseArrays

include("geodesiccomplex.jl")
include("equilateraliterator.jl")
include("persistence.jl")
include("plotting.jl")

export persistence, GeodesicComplex

end
