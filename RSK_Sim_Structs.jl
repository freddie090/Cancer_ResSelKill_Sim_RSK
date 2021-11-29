
# Cancer Resistance Selective Killing Simulation
# RSK.
# Freddie Whiting - 2021

# Mutable Structures

################################################################################


# Cancer Cell - cell_ID, barcode, birth and death rates and resistance score.

mutable struct CancerCell
    cell_ID::Int64 # Cell lineage ID.
    muts::Array{Float64} # Vector of cell's mutations.
    b::Float64 # birth rate.
    d::Float64 # death rate.
    R::Float64 # Resistant phenotype (binary).
    Ndiv::Int64 # Number of divisions a given cell lineage has experienced.
end

# Output of grow_cells.

mutable struct Grow_Out
    cells::Array{CancerCell}
    Nvec::Array{Int64}
    tvec::Array{Float64}
    Rvec::Array{Int64}
    fin_t::Float64
end
