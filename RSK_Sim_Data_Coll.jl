
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

    # Extract all mutations.
    all_muts = vcat(map(x -> x.muts, cell_vec)...)
    # Convert to count dataframe
    mut_df = DataFrame(mut = all_muts)
    mut_df = combine(groupby(mut_df, :mut), nrow)
    colnames = [:mut_ID, :mut_count]
    rename!(mut_df, colnames)
    # Add a relative frequency column
    mut_df[!, :mut_rf] = mut_df[!, :mut_count]./length(cell_vec)

    return mut_df

end
