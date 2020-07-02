using CSV
include("./stochastic_model.jl")
using .stochastic_model
using DataFrames

# parameters
n_cells = 2000
# cell can be alive or dead
# create cell array
time_arr = range(0, 5, step = 0.05)

# first entry is cell cell_state
# second entry is jump time (last state transition)
# third entry is cumulative
function run_sim(n_sim, n_cells, time_arr)
    res_arr = [stoc_model(n_cells, time_arr) for i = 1:n_sim]
    res_arr = vcat(res_arr...)
    return res_arr
end

res_arr = run_sim(100, n_cells, time_arr)
CSV.write("../output/teststoring2.csv", res_arr)
