
# Cancer Resistance Selective Killing Simulation
# RSK.
# Freddie Whiting - 2021

# Data Collection

################################################################################

# Replace missings in a DataFrame.

function rep_missings(df::DataFrame)

    for col in names(df)
        df[ismissing.(df[:,col]), col] .= 0.0
    end

    return df

end


# Full join an array of dataframes by colname and return, removing missings.

function join_dfs(dfs::Array{DataFrame}, colname::String)

    length(dfs) >= 2 || error("Provide >1 DataFrame in the dfs array.")

    full_df = dfs[1]
    for i in 2:length(dfs)
        full_df = outerjoin(full_df, dfs[i], on = Symbol(colname))
    end
    full_df = rep_missings(full_df)

    return full_df

end


# Extract the frequency of mutations from a Grow_Out structure.

function extract_mut_freqs(cell_vec::Array{CancerCell})

    # Extract all mutations - rows = mutation ID, cols = freq.
    all_muts = hcat(map(x -> x.muts, cell_vec)...)
    # Convert into freq vector
    mf_vec = sum(all_muts, dims = 2)
    # Only keep the positive values - others are extinct mutations.
    mf_vec = mf_vec[mf_vec .> 0]
    # Convert into a freq dataframe
    mut_df = DataFrame(mut_ID = collect(1:length(mf_vec)),
                       mut_count = mf_vec,
                       mut_rf = mf_vec/length(cell_vec))

    return mut_df

end
