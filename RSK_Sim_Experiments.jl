
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


# Run a grow-kill-grow-VAF experiment with the given parameters and save the
# output.

function RSK_Run_Exp(nmut::Int64, mutmax::Int64, b::Float64, d::Float64,
    tmax_1::Float64, tmax_2::Float64, Nmax_1::Int64, Nmax_2::Int64, mu::Float64,
    r::Float64, tol_dis::Float64, ploidy::Int64, cellularity::Float64,
    detec_lim::Float64, depth::Int64, Nsim::Int64)

    cd("Outputs")

    # Check for dir, and create if not
    dir_name = string("RSK_Exp_out_nmut-", nmut, "_b-", round(b, digits=3),
    "_d-", round(d, digits = 3), "_tmax_1-", round(tmax_1, digits = 10),
    "_tmax_2-", round(tmax_2, digits = 10), "_Nmax_1-", Nmax_1,
    "_Nmax_2-", Nmax_2, "_mu-", round(mu, digits = 10), "_r-",
    round(r, digits = 10), "_tol_dis-", round(tol_dis, digits = 10),
    "_ploidy-", ploidy, "_cellu-", round(cellularity, digits = 10),
    "_dtclm-", round(detec_lim, digits = 10), "_depth-", depth)


    if isdir(dir_name) == true
        cd(dir_name)
    else
        mkdir(dir_name)
        cd(dir_name)
    end

    NameSim = string("Sim_", Nsim)
    # Make a sub-dir for the different sim iterations
    if isdir(NameSim) == true
        cd(NameSim)
    else
        mkdir(NameSim)
        cd(NameSim)
    end

    # Run the grow-kill-grow-VAF function with the given parameters.
    VAF_df = grow_kill_grow_VAF(nmut, mutmax, b, d,
        tmax_1, tmax_2, Nmax_1, Nmax_2, mu,
        r, tol_dis; ploidy=ploidy, cellularity=cellularity,
        detec_lim=detec_lim, depth=depth)

    @rput VAF_df
    R"""
    write.csv(VAF_df, file = "VAF_df.csv", row.names = F)
    """

    cd("../../")

end


################################################################################
