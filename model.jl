using Flux

struct GraphEmb
    emb::Matrix
    layer::Flux.Dense
end

function GraphEmb(;
        graph_dim::Int,
        emb_dim::Int
        )

    emb = randn(graph_dim, emb_dim)
    layer = Flux.Dense(emb_dim, emb_dim, Flux.relu)
    return GraphEmb(emb, layer)
end

Flux.@functor GraphEmb

function (m::GraphEmb)(x)
    ## if this is a single sample we add a dimention to the end
    if ndims(x) == 2
        x = reshape(x, (size(x)..., 1))
    end
    ## embed the sequence of coords
    x = Flux.batched_mul(x, m.emb)
    ## sum along sequence length
    x = sum(x, dims=1)
    ## drop singleton dim
    x = dropdims(x, dims=1)
    ## pass the summed embedded vector through layer
    x = m.layer(x)
    return x
end

function GraphNet(;
        fc_dim::Int,
        h_dim::Int
        )

    model = Chain(
        Parallel(
            +,
            GraphEmb(
                graph_dim=2, ## (x, y) grid coordinates
                emb_dim=h_dim
            ),
            Dense(fc_dim, h_dim)
        ),
        Dense(h_dim,h_dim)
    )
    return model
end

model = GraphNet(fc_dim=10, h_dim=32)
