
# Cancer Resistance Selective Killing Simulation
# RSK.
# Freddie Whiting - 2021

# Functions

################################################################################

###########################
# Description
###########################

# Description of this version of the Resistance Selective Killing (RSK)
# simtulation.

# Seed a single cell with given birth (b) and death (d) rates. Assign the cell
# an initial number of mutations - these are the clonal mutations found in the
# first cell at the time of malignant transformation.

# 1. Growing cell step:

# Grow the cell until either t = tmax or N = Nmax. Allow mutations to accrue
# with probability mu per-division. Each time a new mutation arises, assign
# it a new mutation identifier.

# 2. Selective killing step:

# Survival is under genetic control:
# Choose a mutation that is found at frequency r. Selectively kill all cells
# that do NOT have this mutation - can therefore think of this as a 'resistance
# mutation' in a drug-treatment setting.
# Survival is not under genetic control:
# In parallel, run a seperate process where we kill all cells with probability
# (1-r), independent of any mutations.
# Can then go on to compare the VAFs of these different scenarios.
# NB that while also want to keep track of the absolute number of cells before
# and after the killing step - assumption is bottleneck size will also have a
# large impact on the VAF.

# 3. Simulate sequencing step:

# Using a ploidy, estimate of cellularity, sequencing depth - convert the
# mutation frequencies into VAFs for quantitative comparisons.

############################
# CopyCell Function.
############################

# Opposed to deepcopying the cell following a birth every time, it is much
# quicker to create a seperate copycell function, and call this within
# the birth-death function.

function copycell(orig_cell::CancerCell)
    new_cell::CancerCell = CancerCell(
    copy(orig_cell.muts),
    copy(orig_cell.b),
    copy(orig_cell.d),
    copy(orig_cell.R),
    copy(orig_cell.Ndiv)
    )
end



############################
# Seed Cells Function.
############################

# Seed a cell with nmut mutations, birth and death rates b and d, respectively.
# Assign the cell with a vector of 0s for each potential mutation position  up
# until mutmax. 1:nmut is populated with 1s. The position in the vector
# corresponds to the mutation identity. NB that in the grow function, if the
# number of unique mutations observed exceeds mutmax, the simulation will end
# and throw an error.

function seed_cell(nmut::Int64, b::Float64, d::Float64)

    # Assign all the fields to the CancerCell accordingly.

    init_cell = CancerCell(collect(1:nmut), b, d, 0, 0)

    return [init_cell]

end


############################
# Growth Function.
############################

# Give an initial cell to the grow cell function. Grow either until t = tmax
# or N = Nmax. Allow mutations to accrue at a rate mu per-division. Record
# the total number of i) cells and ii) mutations every (tmax/t_frac).

function grow_cells(init_cells::Array{CancerCell}, tmax::Float64,
    Nmax::Int64, mu::Float64; t_frac=0.050::Float64)

    out_cells = deepcopy(init_cells)

    0 <= mu || error("mu must be greater than 0.")
    0 <= t_frac <= 1.0 || error("t_frac must be between 0 and 1.")

    bmax = maximum(map(x -> x.b, out_cells))
    dmax = maximum(map(x -> x.d, out_cells))

    length(unique(map(x -> x.b, out_cells))) == 1 || error("For this simulation, cells should have a uniform base birth and death rate.")
    length(unique(map(x -> x.d, out_cells))) == 1 || error("For this simulation, cells should have a uniform base birth and death rate.")

    # Save the population's net growth rate, λ (lam).
    lam = bmax - dmax
    # Save the sum of the rates to normalise the time units.
    bdmax = bmax + dmax

    # Create vector to hold time, and initiate at 0.
    t = 0.0
    tvec = Float64[t]
    # Record the population size (Nt) the total number of mutation mutations
    # (Numuts) every tmax*t_frac.
    t_rec_change = tmax*t_frac
    t_rec = t_rec_change

    # Have a seperate counter that keeps track of the number of live cells,
    # as we can't use length(out_cells) - includes the dead cells until the
    # end.
    Nt = length(out_cells)
    # Have a vector to save the number of live cells every t_rec.
    Nvec = Int64[Nt]
    # Also have a vector to save the number of unique mutations (Numuts).
    Numut = length(unique(vcat(map(x -> x.muts, out_cells)...)))
    # And the most recent unique mutation identifier (not the same as unique
    # mutations, as some mutant lineages could be lost to drift).
    if Numut == 0
        MRmut = 0
    else
        MRmut = maximum(unique(vcat(map(x -> x.muts, out_cells)...)))
    end

    Numutvec = Int64[Numut]

    # Opposed to sampling directly from out_cells, have a 'samp_vec'. Update
    # this following a birth or death accordingly:
    # if a birth, +1 to the length of the vector.
    # if a death, change this index in the sampling vec to 0, and add the index
    # to the 'kill_vec'.
    # Keep sampling from samp_vec until ran_cell_pos > 0
    # This way, don't have to use deleteat! (which is expensive) until the very
    # end.
    samp_vec = collect(1:length(out_cells))
    kill_vec = Array{Int64}(undef, 0)

    # Grow cells until t = tmax

    while t <= tmax

        ran_samp_pos = rand(1:length(samp_vec))
        ran_cell_pos = samp_vec[ran_samp_pos]
        # Keep sampling from samp_vec until hit a +ve number.
        if ran_cell_pos == 0
            while ran_cell_pos == 0
                ran_samp_pos = rand(1:length(samp_vec))
                ran_cell_pos = samp_vec[ran_samp_pos]
            end
        end

        ran = rand(Uniform(0, bdmax))

        # Update the time pre-event, to prevent including events that will
        # have occured post-tmax.

        dt = -1/(bdmax * Nt) * log(rand())
        t += dt

        # If the time-change has crossed t_rec, then update the tvec and
        # Nvec tracking vectors, and update t_rec by t_rec_change.
        if t >= t_rec
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Numutvec, Numut)
            t_rec += t_rec_change
        end

        # If the time-change has crossed tmax, then there wasn't time for
        # the birth or event to happen, t is capped at tmax, and the
        # simulation ends.
        if t > tmax
            t = tmax
            break
        end
        # If N is now >= Nmax, also end the simulation, whilst updating the
        # tracking vectors.
        if Nt >= Nmax
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Numutvec, Numut)
            break
        end

        # Continue with the birth and death steps.
        cell_b = out_cells[ran_cell_pos].b
        cell_d = out_cells[ran_cell_pos].d

        cell_b > 0.0 || error("A cell's birth rate (b) should be positive.")
        cell_d >= 0.0 || error("A cell's death rate (d) should not be negative.")

        # Birth step
        if ran < cell_b

            # Update the cells number of divisions.
            out_cells[ran_cell_pos].Ndiv += 1
            ran_cell = copycell(out_cells[ran_cell_pos])

            # Cells gain a mutation with rate mu.
            # NB Here, both daughter cells inherit the mutation.
            # 'Grey -> Black + Black'

            # Draw number of mutations from a Poisson distribution.
            n_mut = rand(Poisson(mu))

            if n_mut > 0
                # Vector of new mutations
                n_mut_vec = collect((MRmut+1):(MRmut+n_mut+1))
                # Add these mutations to the mother and daughter cells.
                append!(ran_cell.muts, n_mut_vec)
                append!(out_cells[ran_cell_pos].muts, n_mut_vec)
                # Increase the number of unique mutations and the most recent
                # unique mutation identifier in the simulation.
                Numut += (n_mut+1)
                MRmut += (n_mut+1)

            end

            # Now add the daughter cell to the vector of cells.
            push!(out_cells, ran_cell)
            # Update samp_vec with the new length of the output cells.
            push!(samp_vec, (length(out_cells)))
            # Update number of cells to normalise dt.
            Nt += 1
        end

        # Death step
        if bmax <= ran < bmax + cell_d
            # Remove this chosen cell's index from the samp_vec.
            samp_vec[ran_samp_pos] = 0
            push!(kill_vec, ran_cell_pos)
            # Update number of cells to normalise dt.
            Nt -= 1
        end

        # Break if no cells remain. Make use of the fact that:
        # samp_vec = N0 + nbirths
        # kill_vec = ndeaths
        # therefore if length(kill_vec) >= length(samp_vec), there are no cells
        # remaining. Also save this outcome to the record vectors.
        if length(kill_vec) >= length(samp_vec)
            push!(Nvec, Nt)
            push!(tvec, t)
            push!(Numutvec, Numut)
            break
        end

    end

    # Now perform all the kill steps using the indexes stored in kill_vec.
    deleteat!(out_cells, sort(kill_vec))

    fin_t = round(t, digits = 4)

    return Grow_Out(out_cells, Nvec, tvec, Numutvec, fin_t)

end


############################
# Simulate VAF Function.
############################

function sim_VAF(mut_df::DataFrame, mut_rf_col::String; ploidy=2::Int64,
    cellularity=1.0::Float64, detec_lim=0.10::Float64, depth=100::Int64)

    0 <= cellularity <= 1.0 || error("cellularity must be between 0.0 and 1.0.")
    0 <= detec_lim <= 1.0 || error("detec_lim must be between 0.0 and 1.0.")

    # Extract the chosen mutation relative frequencies.
    mut_rfs = mut_df[! , Symbol(mut_rf_col)]
    # Adjust for ploidy
    mut_rfs = mut_rfs./ploidy
    # Adjust for cellularity
    mut_rfs = mut_rfs.*cellularity
    # Set those < detection limit to 0.0
    mut_rfs[mut_rfs .< detec_lim] .= 0.0
    # Sample depths per allele using a binomial where rate = depth.
    all_dep = rand(Poisson(depth), length(mut_rfs))
    # Now sample alleles where n = sampled all_dep, and p = mut_rfs.
    samp_all = map((n, p) -> rand(Binomial(n, p)), all_dep, mut_rfs)
    # And depth is these counts normalsed by the total depth.
    VAF = samp_all ./ all_dep
    # Return the vector of simulated VAFs

    return VAF

end



############################
# Selective Killing Function
############################

# Choose a 'resistant mutation' frequency, r, and find a mutation with this
# given frequency (+/- some tolerance parameter). Then either:
#   - genetic bottleneck...
#  i) kill all cells without this mutation (killing Nt*(1 - r) cells)
#   - non-genetic bottleneck...
# ii) kill all cells with probability (1 - r)
# Can then compare the VAF spectrum given these two selective killing scenarios.

function selectively_kill(cell_vec::Array{CancerCell}, r::Float64,
    tol_dis::Float64; gen_bottle=true::Bool)

    0.0 <= r <= 1.0 || error("r must be between 0.0 and 1.0")

    if gen_bottle == true
        # Genetic bottleneck killing:
        # Calculate the mutation frequency table.
        mut_df = extract_mut_freqs(cell_vec)
        # Extract mutation relative frequencies
        mut_rf = mut_df[!,:mut_rf]
        # Look for mutation locations with frequency r +/- detec lim
        chosen_muts = abs.(mut_rf .- r) .< tol_dis
        # Throw error if there weren't any
        sum(chosen_muts) > 0 || error("There were no mutations with frequency 'r +/- tol_dis'. Try a less stringent tolerance.")
        # Retrieve locations of chosen mutations.
        chosen_mut_locs = findall(x -> x.==1, chosen_muts)
        # Randomly choose one.
        ran_mut_row = sample(chosen_mut_locs)
        # Therefore exract the mut_ID from the VAF dataframe
        ran_mut = mut_df[ran_mut_row, :mut_ID]
        # Get the positions in the cell vector of all that have this mutation.
        cell_mut_locs = findall(x -> sum(x.muts .== ran_mut).==1, cell_vec)
        new_cell_vec = cell_vec[cell_mut_locs]
    else
        # Non-genetic bottleneck killing:
        # Randomly select r cells to keep (/1-r cells are 'killed')
        surv_n = Int64(ceil(length(cell_vec)*r))
        new_cell_vec = cell_vec[sample(1:length(cell_vec), surv_n)]

    end

    return new_cell_vec

end


#############################
# Grow-Kill-Grow Function
#############################

# Grow cells, selectively kill, grow again, then create the mutation frequency
#  dataframe.
# The aim is to simulate something akin to: an initial expansion stage
# following transformation, a selective bottleneck - for example
# drug-treatment, and then a subsequent expansion stage (where mutations
# continue to accumulate). Finally, the mutation frequency distribution
# is returned.
# Always repeat the selective killing step for both a genetic and non-genetic
# bottleneckon the same counts.

# To aid in comparisons, allow multiple values of r to be run on the
# same original cell output... this way can exclude differences due to
# stochasticity in the simulations when making comparisons.

function grow_kill_grow(nmut::Int64, b::Float64, d::Float64,
    tmax_1::Float64, tmax_2::Float64, Nmax_1::Int64, Nmax_2::Int64, mu::Float64,
    rs::Array{Float64, 1}, tol_diss::Array{Float64, 1})

    # First cell:
    init_cell = seed_cell(nmut, b, d)

    # Grow for first interval - repeat if lost to extinction.
    out_1 = grow_cells(init_cell, tmax_1, Nmax_1, mu)

    while length(out_1.cells) == 0
        out_1 = grow_cells(init_cell, tmax_1, Nmax_1, mu)
    end

    # Repeat the following for each value of depth and r, saving the original
    # mutation frequency distribution and selective killing mutation frequency distribution
    # dataframes each time.
    mf_dfs = Array{DataFrame}(undef, 0)

    # Check length of rs = length of tol_diss
    length(rs) == length(tol_diss) || error("The vectors 'rs' and 'tol_diss' should be the same length.")

    for i in 1:length(rs)

        # Save this as an original mut_freq_df.
        orig_mf_df = extract_mut_freqs(out_1.cells)

        # Selectively kill:
        # Genetic bottleneck -
        sel_kill_out_a = selectively_kill(out_1.cells, rs[i], tol_diss[i], gen_bottle=true)
        # Non-genetic bottleneck -
        sel_kill_out_b = selectively_kill(out_1.cells, rs[i], tol_diss[i], gen_bottle=false)

        # Repeat the growth period - again, repeat if lost to extinction.
        # Do this for both genetic and non-genetic bottlenecks.

        out_2a = grow_cells(sel_kill_out_a, tmax_2, Nmax_2, mu)

        while length(out_2a.cells) == 0
            out_2a = grow_cells(sel_kill_out_a, tmax_2, Nmax_2, mu)
        end

        out_2b = grow_cells(sel_kill_out_b, tmax_2, Nmax_2, mu)

        while length(out_2b.cells) == 0
            out_2b = grow_cells(init_cells, tmax_2, Nmax_2, mu)
        end

        # Turn these final outputs into mutation frequency dataframes.

        mf_df_a = extract_mut_freqs(out_2a.cells)

        mf_df_b = extract_mut_freqs(out_2b.cells)

        # Merge all three mutation frequency dataframes:
        mf_df = outerjoin(orig_mf_df, mf_df_a, on = :mut_ID, makeunique=true)
        mf_df = outerjoin(mf_df, mf_df_b, on = :mut_ID, makeunique=true)
        # Replace NAs
        mf_df = rep_missings(mf_df)

        # Rename columns
        colnames = [:mut_ID, :mut_count_pre, :mut_rf_pre,
                             :mut_count_gb, :mut_rf_gb,
                             :mut_count_ngb, :mut_rf_ngb]
        rename!(mf_df, colnames)

        # Add column saving this simulations given r, and then save to
        # the output arrays.

        mf_df[!, :r] = repeat([rs[i]], nrow(mf_df))
        mf_df[!, :tol_dis] = repeat([tol_diss[i]], nrow(mf_df))

        push!(mf_dfs, mf_df)

    end

    # Merge for each r and depth
    mf_df = vcat(mf_dfs...)

    # Return final dataframes
    return mf_df

end
