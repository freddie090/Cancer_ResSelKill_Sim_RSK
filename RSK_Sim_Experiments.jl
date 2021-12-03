
# Cancer Barcode Simulation
# Freddie Whiting - 2021

# Experiment Functions

# A collection of functions that run pre-made experimental set-ups.

# To distinguish from other functions, all those that run an experiment should
# begin with 'Run_Exp...'.

################################################################################

# Import Threads functions for multi-threading.
import Base.Threads.@threads
import Base.Threads.@spawn


# Experimental set-up.

# 4 Replicates per treatment (Control vs Drug-Treatment).
# If Passage number > 1, can provide the output of a previous Passage as
# an 'Experiment_Input' object to continue with a N_seed randomly sampled
# cells from the previous flask.

function Run_Exp(N::Int64, b::Float64, d::Float64, p::Float64, mu::Float64,
    sig::Float64, del::Float64, R_real::String,
    n_pulse::Int64, Nmax::Int64, N_seed::Int64,
    t_CO::Float64, t_DT::Float64, Passage::Int64,
    insta_kill::Bool, lim_probs::Bool;
    al=0.0::Float64, psi=0.0::Float64,
    Exp_Input::Experiment_Input=Experiment_Input())

    if Passage == 1

        # First expand and split the cells for t = 6.0.
        #exp_split_cells = expand_split_cells(N, b, d, 6.0, N_seed, 4, p, mu, sig, del, psi=psi, al=al, R_real=R_real, use_lim_probs=lim_probs)
        # First expand and split the cells for t = 7.0: for QR_HCTbc.
        # exp_split_cells = expand_split_cells(N, b, d, 7.0, N_seed, 4, p, mu, sig, del, psi=psi, al=al, R_real=R_real, use_lim_probs=lim_probs)
        # First expand and split the cells for t = 9.0: for QR_SW6bc.
        exp_split_cells = expand_split_cells(N, b, d, 9.0, N_seed, 4, p, mu, sig, del, psi=psi, al=al, R_real=R_real, use_lim_probs=lim_probs)


        Exp_Input.CO_inputs = exp_split_cells.CO_flasks
        Exp_Input.DT_inputs = exp_split_cells.DT_flasks
    end

    CO_Outputs = Array{Grow_Kill_Rec_Out}(undef, 0)
    DT_Outputs = Array{Grow_Kill_Rec_Out}(undef, 0)

    c1 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,1], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 1, "P", Passage)))
    c2 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,2], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 2, "P", Passage)))
    c3 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,3], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 3, "P", Passage)))
    c4 = @spawn push!(CO_Outputs, grow_kill_rec_cells(Exp_Input.CO_inputs[:,4], t_CO, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, rep_name = string("CO", 4, "P", Passage)))
    d1 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,1], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 1, "P", Passage)))
    d2 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,2], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 2, "P", Passage)))
    d3 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,3], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 3, "P", Passage)))
    d4 = @spawn push!(DT_Outputs, grow_kill_rec_cells(Exp_Input.DT_inputs[:,4], t_DT, mu, sig, del, n_pulse, Nmax, R_real=R_real, psi=psi, al=al, drug_kill=true, insta_kill=insta_kill, rep_name = string("DT", 4, "P", Passage)))

    wait(c1); wait(c2); wait(c3); wait(c4); wait(d1); wait(d2); wait(d3); wait(d4);

    return Experiment_Output(CO_Outputs, DT_Outputs)

end


# Run the experiment and save the output.

function Run_Exp_save_output(N::Int64, b::Float64, d::Float64, p::Float64,
    mu::Float64, sig::Float64, del::Float64,
    R_real::String, n_pulse::Int64, Nmax::Int64,
    N_seed::Int64, t_CO::Float64, t_DT::Float64, Nsim::Int64, Passage::Int64,
    insta_kill::Bool, lim_probs::Bool;
    psi=0.0::Float64, al=0.0::Float64)

    cd("Outputs")

    # Save insta_kill and lim_probs as either 0/1.

    if insta_kill == true
        insta_kill_state = "1"
    elseif insta_kill == false
        insta_kill_state = "0"
    end
    if lim_probs == true
        lim_probs_state = "1"
    elseif lim_probs == false
        lim_probs_state = "0"
    end


    # Check for dir, and create if not
    dir_name = string("CBC_Exp_out_N-", N, "_b-", round(b, digits=3),
    "_d-", round(d, digits = 3), "_p-", round(p, digits = 10),
    "_mu-", round(mu, digits = 10), "_sig-", round(sig, digits = 10),
    "_del-", round(del, digits = 10),
    "_psi-", round(psi, digits = 10), "_al-", round(al, digits = 10),
    "_R_real-", R_real,
    "_n_pulse-", n_pulse, "_Nmax-", Nmax,
    "_N_seed-", N_seed, "_t_CO-", t_CO, "_t_DT-", t_DT,
    "_ik-", insta_kill_state, "_lp-", lim_probs_state)

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

    # Call outside the for loop first for scoping.
    exp_in = Experiment_Input()

    for P in 1:Passage
        if P == 1
            exp_out = Run_Exp(N, b, d, p, mu, sig, del, R_real, n_pulse, Nmax,
            N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
            psi=psi, al=al)
        else

            exp_out = Run_Exp(N, b, d, p, mu, sig, del, R_real, n_pulse, Nmax,
            N_seed, t_CO, t_DT, P, insta_kill, lim_probs,
            psi=psi, al=al,
            Exp_Input=exp_in)
        end

        # If the current Passage has a replicate that has no cells
        # remaining, kill this sim.
        if sum(map(x -> length(x.cells), exp_out.CO_outputs) .== 0) > 0
            println("A control replicate has 0 cells remaining. Ending sim.")
            break
        end
        if sum(map(x -> length(x.cells), exp_out.DT_outputs) .== 0) > 0
            println("A drug-treatment replicate has 0 cells remaining. Ending sim.")
            break
        end
        # Save both the total number of resistant cells (R > 0.0) and the mean
        # resistant phenotype (mean(R)) as well as the counts per barcode
        # lineage.
        bc_R_df = all_bc_counts_tot_and_mean_R(exp_out)
        #NRE_df = all_NRE_by_ts(exp_out)
        N_df = all_N_by_ts(exp_out)

        #@rput bc_R_df ; @rput NRE_df ; @rput P
        @rput bc_R_df; @rput N_df; @rput P
        R"""
        write.csv(bc_R_df, file=paste0("bc_counts_tot_mean_R_P", P, ".csv"), row.names = F)
        #write.csv(NRE_df, file=paste0("NRE_by_t_P", P, ".csv"), row.names = F)
        write.csv(N_df, file=paste0("N_by_t_P", P, ".csv"), row.names = F)
        """

        # Might have had enough to save the barcode distributions and N vector,
        # but not enough cells to draw the next Passage's flasks from. Check
        # for that now.
        if sum(map(x -> length(x.cells), exp_out.CO_outputs) .< N_seed) > 0
            println("Not enough cells in the control replicates to seed the next Passages flask. Ending sim.")
            break
        end
        if sum(map(x -> length(x.cells), exp_out.DT_outputs) .< N_seed) > 0
            println("Not enough cells in the drug-treatment replicates to seed the next Passages flask. Ending sim.")
            break
        end


        exp_in = Experiment_Input()
        exp_in.CO_inputs=Array{CancerCell,2}(undef, N_seed, 4)
        exp_in.DT_inputs=Array{CancerCell,2}(undef, N_seed, 4)
        for i in 1:4
            exp_in.CO_inputs[:,i] = sample(exp_out.CO_outputs[i].cells, N_seed, replace=false)
            exp_in.DT_inputs[:,i] = sample(exp_out.DT_outputs[i].cells, N_seed, replace=false)
        end

    end

    cd("../../../")

end


################################################################################
