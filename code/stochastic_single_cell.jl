# use this script to test rate/poisson/resposne time stuff for a single cell

using Distributions
using DataFrames
using CSV

# for each time point loop over each cell

function naive_diff(cell_arr, i, time, alpha, beta, th1_idx, cell_state, last_state_change,
prob_state_change, prob_death, prob_div, prob_div0, last_death, last_div)
# note that naive cell does not need to take jump instance into account because last_state change = 0
    if 1-cdf(Gamma(alpha,1/beta), time) < cell_arr[i,prob_state_change]
        # assign new cell state
        cell_arr[i,cell_state] = th1_idx

        # assign jump time
        cell_arr[i,last_state_change] = time
        cell_arr[i,last_div] = time
        cell_arr[i, last_death] = time
        # add this only because only effecotr cells die
        cell_arr[i, prob_state_change] = rand(Float64)
    end
end


function cell_prolif(cell_arr, i, time, t_new, r_div, dead_idx, cell_state, last_state_change,
prob_state_change, prob_death, prob_div, prob_div0, last_death, last_div)

    if cell_arr[i, prob_div0] > cell_arr[i, prob_div]

        create_cell(cell_arr, i, time, dead_idx, cell_state, last_state_change,
        prob_state_change, prob_death, prob_div, prob_div0, last_death, last_div)

        # for mother cell update division time and assign new division probability
        cell_arr[i, prob_div] = rand(Float64)
        cell_arr[i, prob_div0] = 0
        cell_arr[i, last_div] = time
    else
        cell_arr[i, prob_div0] = cell_arr[i, prob_div0] +
        cdf(Exponential(1/r_div),t_new-cell_arr[i,last_div])-
        cdf(Exponential(1/r_div),time-cell_arr[i,last_div])
    end
end

function create_cell(cell_arr, i, time, dead_idx, cell_state, last_state_change,
prob_state_change, prob_death, prob_div, prob_div0, last_death, last_div)
    # find a dead cell and make it alive
    # from all dead cells start with the first index then increase
    dead_cells = findall(cell_arr[:, cell_state] .== dead_idx)
    k = dead_cells[1]
    cell_arr[k,cell_state] = cell_arr[i, cell_state]
    # for daughter cell add new div and death probabilities
    # also assign new last death/div/state change time
    cell_arr[k, prob_div] = rand(Float64)
    cell_arr[k, prob_death] = rand(Float64)
    # no one cares about state change at this point but anyway
    cell_arr[k,last_state_change] = time
    cell_arr[k,last_div] = time
    cell_arr[k,last_death] = time
    cell_arr[k, prob_div0] = 0
end


function cell_death(cell_arr, i, time, r_death, dead_idx, cell_state, last_state_change,
prob_state_change, prob_death, prob_div, prob_div0, last_death, last_div)

    if (1-(cdf(Exponential(1), time-cell_arr[i,last_death]))) < cell_arr[i,prob_death]
        cell_arr[i,cell_state] = dead_idx
        cell_arr[i, last_death] = time
    end
end

function stoc_model(n_cells, time_arr)
    ########################################################## assign parameters
    alpha = 2
    beta = 2
    r_div = 1.0
    r_death = 1.0
    ############################################################# assign indices
    naive_idx = 0
    th1_idx = 1
    dead_idx = 8
    ######################################################### assign cell states
    cell_state = 1
    last_state_change = 2
    prob_state_change = 3
    prob_death = 4
    prob_div = 5
    prob_div0 = 6
    last_death = 7
    last_div = 8

    cell_params = [cell_state, last_state_change, prob_state_change, prob_death,
    prob_div, prob_div0, last_death, last_div]

    ######################################################### create cell array
    n_states = length(cell_params)
    cell_arr = zeros((n_cells, n_states))
    # assign random numbers for naive prec. transition
    cell_arr[:,prob_state_change] = rand(n_cells)
    # assign random numbers for death event of cells
    cell_arr[:,prob_death] = rand(n_cells)
    cell_arr[:,prob_div] = rand(n_cells)

    n_dead_cells = 1000
    dead_cell_arr = zeros((n_dead_cells, n_states))
    dead_cell_arr[:, cell_state] .= dead_idx
    cell_arr = [cell_arr; dead_cell_arr]

    ################################################ create arr for cell numbers
    n_naive = zeros(length(time_arr))
    n_th1 = zeros(length(time_arr))

    for t = 1:(length(time_arr)-1)
        time = time_arr[t]
        t_new = time_arr[t+1]
        # get indices of dead cells
        # use 8 as idx for dead cells
        #alive_cells = findall(cell_arr[:,cell_state] .!= dead_idx)

        for i = 1:(n_cells+n_dead_cells)
            # note that the order of death vs prolif matters here
            # if they cells die first they cannot proliferate, reverse is not true
            # still not sure if I should use update state at t+1 so that all effects

            # differentiation
            if cell_arr[i,cell_state] == naive_idx
                naive_diff(cell_arr, i, time, alpha, beta, th1_idx, cell_params...)
            end

            # proliferation
            if cell_arr[i,cell_state] in [th1_idx]
                if cell_arr[i, prob_div0] > cell_arr[i, prob_div]

                    dead_cells = findall(cell_arr[:, cell_state] .== dead_idx)
                    k = dead_cells[1]
                    cell_arr[k,cell_state] = cell_arr[i, cell_state]
                    # for daughter cell add new div and death probabilities
                    # also assign new last death/div/state change time
                    cell_arr[k, prob_div] = rand(Float64)
                    cell_arr[k, prob_death] = rand(Float64)
                    # no one cares about state change at this point but anyway
                    cell_arr[k,last_state_change] = time
                    cell_arr[k,last_div] = time
                    cell_arr[k,last_death] = time
                    cell_arr[k, prob_div0] = 0

                    # for mother cell update division time and assign new division probability
                    cell_arr[i, prob_div] = rand(Float64)
                    cell_arr[i, prob_div0] = 0
                    cell_arr[i, last_div] = time
                else
                    cell_arr[i, prob_div0] = cell_arr[i, prob_div0] +
                    cdf(Exponential(1/r_div),t_new-cell_arr[i,last_div])-
                    cdf(Exponential(1/r_div),time-cell_arr[i,last_div])
                end
            end

            # death
            if cell_arr[i,cell_state] in [th1_idx]
                if (1-(cdf(Exponential(1/r_death), time-cell_arr[i,last_death]))) < cell_arr[i,prob_death]
                    cell_arr[i,cell_state] = dead_idx
                    cell_arr[i, last_death] = time
                end
            end
        end

        # get cell numbers
        n_th1[t] = sum(cell_arr[:,cell_state] .== th1_idx) #println(t)
    end

    df = DataFrame(time = time_arr, Th1 = n_th1)
    return df
end

#################################################### change below if I want module
#export stoc_model
#end
# check how update rules for cells will apply


# parameters
n_cells = 500
# cell can be alive or dead
# create cell array
time_arr = range(0, 20, step = 0.001)

# first entry is cell cell_state
# second entry is jump time (last state transition)
# third entry is cumulative
function run_sim(n_sim, n_cells, time_arr)
    res_arr = [stoc_model(n_cells, time_arr) for i = 1:n_sim]
    res_arr = vcat(res_arr...)
    return res_arr
end

df = run_sim(50, n_cells, time_arr)

#using StatsPlots
#@df df plot(:time, [:Th1 :Tfh], colour = [:red :blue])

CSV.write("Onedrive/projects/2020/tcell_stochastic/output/single_cell.csv", df)
