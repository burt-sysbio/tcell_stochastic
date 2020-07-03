##################################### add if I want module module stochastic_model

using Distributions
using DataFrames
using CSV

function draw_fate(p1=0.5, p2=0.2, p3=0.3)
    # draw random number
    n = rand(Float64)
    if 0 <= n < p1
        fate = 1
    elseif p1 <= n < p1+p2
        fate = 2
    else
        fate = 3
    end
    return fate
end


function pos_fb(c, EC50)
    return c/(c+EC50)
end

function get_myc(time, deg_myc)
    return exp(-deg_myc*time)
end
# for each time point loop over each cell
function stoc_model(n_cells, time_arr)

    alpha = 2
    beta = 2
    r_div_base = 1
    EC50_myc = 0.5
    deg_myc = 0.5
    # assign indices
    naive_idx = 0
    th1_idx = 1
    tfh_idx = 2
    tr1_idx = 3
    dead_idx = 8

    cell_state = 1
    last_state_change = 2
    prob_state_change = 3
    prob_death = 4
    prob_div = 5
    prob_div0 = 6
    last_death = 7
    last_div = 8

    n_states = 8
    cell_arr = zeros((n_cells, n_states))
    # assign random numbers for naive prec. transition
    cell_arr[:,prob_state_change] = rand(n_cells)
    # assign random numbers for death event of cells
    cell_arr[:,prob_death] = rand(n_cells)
    cell_arr[:,prob_div] = rand(n_cells)

    n_dead_cells = 2000
    dead_cell_arr = zeros((n_dead_cells, n_states))
    dead_cell_arr[:, cell_state] .= dead_idx
    cell_arr = [cell_arr; dead_cell_arr]

    n_naive = zeros(length(time_arr))
    n_th1 = zeros(length(time_arr))
    n_tfh = zeros(length(time_arr))
    n_tr1 = zeros(length(time_arr))

    for t = 1:(length(time_arr)-1)
        time = time_arr[t]
        t_new = time_arr[t+1]
        myc = get_myc(time, deg_myc)
        r_div = r_div_base*pos_fb(myc, EC50_myc)
        # get indices of dead cells
        # use 8 as idx for dead cells
        alive_cells = findall(cell_arr[:,cell_state] .!= dead_idx)

        for i in alive_cells
            # not taking jump instance into account right now cuz naive cells
            # use so that I can increase it
            if cell_arr[i,cell_state] == naive_idx
                if 1-cdf(Gamma(alpha,1/beta), time) < cell_arr[i,prob_state_change]
                    # assign new cell state
                    cell_arr[i,cell_state] = draw_fate()
                    # assign jump time
                    cell_arr[i,last_state_change] = time
                    cell_arr[i,last_div] = time
                    # add this only because only effecotr cells die
                    cell_arr[i,last_death] = time
                end
            end

            # for each effector cell type
            ###########################
            #############  note that the order of death vs prolif matter here
            # if they cells die first they cannot proliferate, reverse is not true
            # if cells proiferate they can still die afterwards
            # still not sure if I should use update state at t+1 so that all effects
            # can take place irrespective of ordering
            if cell_arr[i,cell_state] == th1_idx
                #if (1-(cdf(Exponential(1/r_div), time-cell_arr[i,last_div]))) < cell_arr[i,prob_div]

                if cell_arr[i, prob_div0] > cell_arr[i, prob_div]    # find a dead cell and make it alive
                    # from all dead cells start with the first index then increase
                    dead_cells = findall(cell_arr[:,cell_state] .== dead_idx)
                    k = dead_cells[1]
                    cell_arr[k,cell_state] = th1_idx
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
            # check if cell state is not naive (0) or dead(8)
            # then check if cell will die
            if cell_arr[i,cell_state] in [th1_idx, tfh_idx, tr1_idx]
                if (1-(cdf(Exponential(1), time-cell_arr[i,last_death]))) < cell_arr[i,prob_death]
                    cell_arr[i,cell_state] = dead_idx
                    # do I need to do anything with dead cells?
                    # when dead cell comes alive I assign everythin new
                    #cell_arr[i,last_death] = time
                    #cell_arr[i,last_state_change] = time
                end
            end

        end

        # get cell numbers
        n_th1[t] = sum(cell_arr[:,cell_state] .== th1_idx)
        n_tfh[t] = sum(cell_arr[:,cell_state] .== tfh_idx)
        n_tr1[t] = sum(cell_arr[:,cell_state] .== tr1_idx)
        n_naive[t] = sum(cell_arr[:,cell_state] .== naive_idx)    #println(t)
    end
    df = DataFrame(time = time_arr, Th1 = n_th1, Tfh = n_tfh, Tr1 = n_tr1)
    return df
end

#################################################### change below if I want module
#export stoc_model
#end
# check how update rules for cells will apply


# parameters
n_cells = 1000
# cell can be alive or dead
# create cell array
time_arr = range(0, 10, step = 0.05)

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

CSV.write("Onedrive/projects/2020/tcell_stochastic/output/teststoring.csv", df)
