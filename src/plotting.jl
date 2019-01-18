using RecipesBase

function getxyz(pts)
    xs = get.(pts, 1, 0.0)
    ys = get.(pts, 2, 0.0)
    zs = get.(pts, 3, 0.0)
    if !isempty(pts) && length(pts[1]) == 3
        xs, ys, zs
    else
        xs, ys
    end
end

@recipe function plot(gc::GeodesicComplex,
                      res = nothing; only_landmarks = true, graph = true)
    # edges
    if graph
        @series begin
            label := "edges"
            edgepoints = Vector{Float64}[]
            for e in edges(landmark_graph(gc))
                s = points(gc, landmarks(gc, src(e)))
                d = points(gc, landmarks(gc, dst(e)))
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
        getxyz(points(gc, landmarks(gc)))
    end
    # others
    if !only_landmarks
        @series begin
            markersize --> 0.5
            seriestype --> :scatter
            label := "others"
            getxyz(points(gc, setdiff(1:n_points(gc), landmarks(gc))))
        end
    end

    if res ≢ nothing
        for (i, c) in enumerate(reverse(res))
            @series begin
                label := "$i, d=$(round(c.diameter, sigdigits=5)), p=$(round(c.perimeter, sigdigits=5))"
                linewidth := 5
                edgepoints = vcat(points(c), [first(points(c))])
                getxyz(edgepoints)
            end
        end
    end
end
