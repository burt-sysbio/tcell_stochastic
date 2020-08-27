include("stochastic_model.jl")
using HDF5

# cell state dict
d_state = Dict([
("state", 1),
("last_change", 2),
("last_death", 3),
("last_div", 4),
("prob_change", 5),
("prob_div0", 6),
("prob_div", 7),
("prob_death", 8)
])

d_fate = Dict([
("Naive", 0),
("Prec", 1),
("Th1", 2),
("Tfh", 3),
("Tr1", 4),
("Dead", 5)
])

d_param = Dict([
("alpha", 2),
("beta", 2.0),
("r_div", 1.0),
("r_death", 1.0),
("prob_Th1", 0.5),
("prob_Tfh", 0.3),
("prob_Tr1", 0.2)
])

"""
run simulation
"""
function run_sim(n_sim, n_cells, n_genes, time_arr, d_param, d_state, d_fate, mode)

    cell_arr = []
    gene_arr = []

    # run stoc_model several times, each time add gene arr and cell arr to list
    Threads.@threads for i = 1:n_sim
        res = stochastic_module.stoc_model(n_cells, n_genes, time_arr, d_param,
        d_state, d_fate, mode)
        cell_df = res[1]
        gene_df = res[2]
        push!(cell_arr, cell_df)
        push!(gene_arr, gene_df)
    end

    # concatenate cell arr to one large list
    res_cells = vcat(cell_arr...)

    return (res_cells, gene_arr)
end

"""
save file as hdf5 format
"""
function save_file(cell_arr, gene_arr, filename)
    # set working directory
    cd(@__DIR__)
    cd("..")
    cd("output")

    sc_data =h5open("scseq_"*filename*".h5","w")
    for i=1:n_sim
        sc_data[string(i)] = gene_arr[i]
    end
    close(sc_data)

    cell_data =h5open("cell_numbers_"*filename*".h5","w")
    cell_data["cell_data"] = cell_arr
    close(cell_data)
end

# set parameters
n_cells = 2000
# cell can be alive or dead
# create cell array
n_genes = 60
n_sim = 10

time_arr = range(0, 5, step = 0.001)

# run simulation with predetermined and with fixed fate
#cell_arr, gene_arr = run_sim(n_sim, n_cells, n_genes, time_arr, d_param,
#d_state, d_fate, "predetermined")

#save_file(cell_arr, gene_arr, "predetermined_fate")

filename = "fixed_probabilities"
cell_arr, gene_arr = run_sim(n_sim, n_cells, n_genes, time_arr, d_param,
d_state, d_fate, filename)

save_file(cell_arr, gene_arr, filename)
