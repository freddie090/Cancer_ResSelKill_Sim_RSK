
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
    # Work backwards, and only discard all of the 0s until the first positive
    # value - these were redundant positions assigned to the mutation array.
    mf_rev = reverse(mf_vec, dims = 1)[:,1]
    pos_rev_vals = findall(x -> x > 0, mf_rev)
    # First value is the reverse index of the first positive value
    last_pos_val = length(mf_vec) - (pos_rev_vals[1] - 1)
    # Only keep the mutations up until this last positive value.
    mf_vec = mf_vec[1:last_pos_val]
    # Convert into a freq dataframe
    mut_df = DataFrame(mut_ID = collect(1:length(mf_vec)),
                       mut_count = mf_vec,
                       mut_rf = mf_vec/length(cell_vec))

    return mut_df

end
