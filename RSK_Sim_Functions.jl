
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
    copy(orig_cell.cell_ID),
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

function seed_cell(nmut::Int64, b::Float64, d::Float64)

    # Assign all the fields to the CancerCell accordingly.

    init_muts = collect(1.0:nmut)

    init_cell = CancerCell(1, init_muts, b, d, 0, 0)

    return init_cell

end


############################
# Growth Function.
############################

# Give an initial cell to the grow cell function. Grow either until t = tmax
# or N = Nmax. Allow mutations to accrue at a rate mu per-division. Record
# the total number of i) cells and ii) mutations every (tmax/t_frac).

function grow_cell(init_cell::CancerCell, tmax::Float64,
    Nmax::Int64, mu::Float64; t_frac=0.050::Float64)

    out_cells = [deepcopy(init_cell)]

    0 <= mu <= 1.0 || error("mu must be between 0 and 1.")
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
    Numuts = Int64(maximum(init_cell.muts))
    mutvec = Int64[Numuts]

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
            push!(mutvec, Numuts)
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
            push!(mutvec, Numuts)
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

            ran_p = rand()

            # Cells gain a mutation with rate mu.
            # NB Here, both daughter cells inherit the mutation.
            # 'Grey -> Black + Black'

            if mu > ran_p
                # Increase the number of unique mutations in the simulation.
                Numuts += 1
                # Add this new mutation to both daughter cells
                push!(ran_cell.muts, Float64(Numuts))
                push!(out_cells[ran_cell_pos].muts, Float64(Numuts))
                # New mutation = new 'sub-clone'
                ran_cell.cell_ID += 1
                out_cells[ran_cell_pos].cell_ID += 1
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
            push!(mutvec, Numuts)
            break
        end

    end

    # Now perform all the kill steps using the indexes stored in kill_vec.
    deleteat!(out_cells, sort(kill_vec))

    fin_t = round(t, digits = 4)

    return Grow_Out(out_cells, Nvec, tvec, mutvec, fin_t)

end




############################
# Simulate VAF Function.
############################

function sim_VAF(cell_vec::Array{CancerCell}; ploidy=2::Int64,
    cellularity=1.0::Float64, detec_lim=0.05::Float64, depth=100::Int64)

    # Total number of cells:
    Nt = length(cell_vec)
    # Now get the mutation dataframe.
    mut_df = extract_mut_freqs(cell_vec)
    # Extract the allele frequencies
    all_freq = mut_df[!, :mut_rf]
    # Adjust for ploidy
    all_freq = all_freq./ploidy
    # Adjust for cellularity
    all_freq = all_freq.*cellularity
    # Set those < detection limit to 0.0
    #all_freq[all_freq .< detec_lim] .= 0.0
    # Sample depths per allele using a binomial where = Nt and p = depth/Nt
    all_dep = rand(Binomial(Nt, depth/Nt), length(all_freq))

    # Now sample alleles where n = sampled all_dep, and p = all_freq.
    samp_all = map((n, p) -> rand(Binomial(n, p)), all_dep, all_freq)
    # And depth is these counts normalsed by the total depth.
    VAF = samp_all ./ all_dep

    # Return these as part of the dataframe.
    VAF_df = mut_df
    VAF_df[!, :VAF] = VAF

    return VAF_df

end
