
# Cancer Resistance Selective Killing Simulation
# RSK.
# Freddie Whiting - 2021

# Experiment Functions

# A collection of functions that run pre-made experimentals.

# To distinguish from other functions, all those that run an experiment should
# begin with 'RSK_Run_Exp...'.

################################################################################

# Import Threads functions for multi-threading.
#import Base.Threads.@threads
#import Base.Threads.@spawn


# Run a grow-kill-grow experiment with the given parameters and save the
# output.
# NB that the VAF simulations are now performed sepearetly.

function RSK_Run_Exp(nmut::Int64, b::Float64, d::Float64,
    tmax_1::Float64, tmax_2::Float64, Nmax_1::Int64, Nmax_2::Int64, mu::Float64,
    rs::Array{Float64, 1}, tol_diss::Array{Float64, 1}, Nsim::Int64)

    cd("Outputs")

    # Check for dir, and create if not
    dir_name = string("RSK_Exp_out_nmut-", nmut,
    "_b-", round(b, digits=3), "_d-", round(d, digits = 3),
    "_tmax_1-", round(tmax_1, digits = 10), "_tmax_2-", round(tmax_2, digits = 10),
    "_Nmax_1-", Nmax_1, "_Nmax_2-", Nmax_2,
    "_mu-", round(mu, digits = 10))


    if isdir(dir_name) == true
        cd(dir_name)
    else
        mkdir(dir_name)
        cd(dir_name)
    end

    # Run the grow-kill-grow-VAF function with the given parameters.
    mut_df = grow_kill_grow(nmut, b, d,
        tmax_1, tmax_2, Nmax_1, Nmax_2, mu,
        rs, tol_diss)

    @rput mut_df; @rput Nsim;
    R"""
    write.csv(mut_df, file = paste0("mut_df_sim_", Nsim, ".csv"),
                row.names = F)
    """

    cd("../../../")

end


################################################################################

# Simulate a VAF with the output RSK_Run_Exp with the given VAF parameters.

function RSK_Run_Exp_Sim_VAF(out_dir::String, ploidy::Int64,
    cellularity::Float64, detec_lim::Float64, depths::Array{Int64, 1})

    # Change to the given output directory.
    cd("Outputs")
    cd(out_dir)

   # Run and save the VAF dataframes for all mutation dataframes.
   curr_files = readdir()
   for i in 1:length(curr_files)
       if occursin(r"mut_df_.*", curr_files[i])

           mut_df = DataFrame(CSV.File(curr_files[i]))

           VAF_dfs = Array{DataFrame}(undef, 0)

           for depth in depths

               VAF_df = deepcopy(mut_df)

               pre_VAFs = sim_VAF(mut_df, "mut_rf_pre", ploidy=ploidy,
               cellularity=cellularity,
               detec_lim = detec_lim, depth = depth)

               gb_VAFs = sim_VAF(mut_df, "mut_rf_gb", ploidy=ploidy,
               cellularity=cellularity,
               detec_lim = detec_lim, depth = depth)

               ngb_VAFs = sim_VAF(mut_df, "mut_rf_ngb", ploidy=ploidy,
               cellularity=cellularity,
               detec_lim = detec_lim, depth = depth)

               length(pre_VAFs) == length(gb_VAFs) == length(ngb_VAFs) || error("The three simulated VAF vectors should be the same length.")

               # Bind these to the mutation dataframe, along with the given
               # simulated VAF parameters.

               VAF_df[!, :pre_VAF] = pre_VAFs
               VAF_df[!, :gb_VAF] = gb_VAFs
               VAF_df[!, :ngb_VAF] = ngb_VAFs

               VAF_df[!, :detec_lim] = repeat([detec_lim], nrow(VAF_df))
               VAF_df[!, :depth] = repeat([depth], nrow(VAF_df))

               push!(VAF_dfs, VAF_df)

           end

           # Join all of the simulated VAF dataframes.
           VAF_df = vcat(VAF_dfs...)
           mut_df_name = curr_files[i]
           # Write the simulated VAF dataframe in R.
           @rput VAF_df; @rput mut_df_name;
           R"""
           VAF_csv_name <- sub("mut", "VAF", mut_df_name)
           write.csv(VAF_df, file = VAF_csv_name, row.names = F)
           """

       end
   end

   cd("../../")

end
