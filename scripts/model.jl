"""
Creates a `Flux` model which is capable of processing an input composed of a sequence and vector.
The sequence is processed with a `LSTM` cell and the last hidden state is retained. The vector is
processed with a `Dense` layer with a relu activation function. The outputs are concatenated.

# Arguments
- `env::DesignEnvironment`: the environment for which to produce the model
- `h_size::Int`: the output size of both `LSTM` and `Dense` layers

# Example
```
model = lstm_block(env, 128)

x = env |> state |> stack |> model
```

"""
function lstm_block(
        env::DesignEnvironment{<: AbstractDesign, <: AbstractObjective, SequenceVectorState},
        h_size::Int)

    seq_input, vec_input = state_space(env)

    ## Determining the input size to the LSTM
    seq_features, seq_length = size(seq_input)
    ## Input size to the Dense layer
    vec_length = size(vec_input)[1]

    ## Defining the pathway for the sequence
    seq_model = Chain(
        LSTM(seq_features, h_size),
        x -> x[:, :, end])
    ## Simple fully connected layer with relu
    vec_model = Dense(vec_length, h_size, relu)

    ## Process the sequence and the vector in Parallel and then concatenate outputs
    model = Parallel(
        vcat,
        seq_model,
        vec_model)

    return model
end