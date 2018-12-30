using RecipesBase

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

@recipe function plot(res::IntrinsicPersistenceResults;
                      only_landmarks = true,
                      dense = false,
                      dead = false)
    gc = res.complex
    # points
    @series begin
        markersize --> 1.0
        seriestype := :scatter
        label := "landmarks"
        getxyz(landmarks(gc))
    end
    # others
    if !only_landmarks
        @series begin
            markersize --> 0.5
            seriestype --> :scatter
            label := "others"

            getxyz(nonlandmarks(gc))
        end
    end
    # cycles
    if !dense
        for (i, c) in enumerate(reverse(res.cycles))
            @series begin
                label := "$i, d=$(round(c.diameter, sigdigits=5)), p=$(round(c.perimeter, sigdigits=5))"
                linewidth := 5
                edgepoints = vcat(points(gc, c.points), [points(gc, c.points[1])])
                getxyz(edgepoints)
            end
        end
    else
        for (i, c) in enumerate(reverse(res.cycles))
            @series begin
                label := "$i, d=$(round(c.diameter, sigdigits=5)), p=$(round(c.perimeter, sigdigits=5))"
                linewidth := 5
                if length(c.dense) ≥ 2
                    pts = vcat(points(gc, c.dense), [points(gc, c.dense[1])])
                end
                getxyz(pts)
            end
        end
    end
end
