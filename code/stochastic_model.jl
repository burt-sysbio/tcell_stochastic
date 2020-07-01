using Distributions
using Plots

function draw_fate(p1=0.2, p2=0.5, p3=0.3)
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


# parameters
n_cells = 500
# cell can be alive or dead
# create cell array
time_arr = range(0, 5, step = 0.05)

# first entry is cell cell_state
# second entry is jump time (last state transition)
# third entry is cumulative
function run_sim(n_sim, n_cells, time_arr)
    res_arr = [stoc_model(n_cells, time_arr) for i = 1:n_sim]
    return res_arr
end
# for each time point loop over each cell
function stoc_model(n_cells, time_arr)
    cell_arr = zeros((n_cells, 4))

    # assign random numbers for naive prec. transition
    cell_arr[:,3] = rand(n_cells)
    # assign random numbers for death event of cells
    cell_arr[:,4] = rand(n_cells)

    n_naive = zeros(length(time_arr))
    n_th1 = zeros(length(time_arr))
    n_tfh = zeros(length(time_arr))
    n_tr1 = zeros(length(time_arr))

    for t = 1:(length(time_arr)-1)
        time = time_arr[t]

        # get indices of dead cells
        # use 8 as idx for dead cells
        dead_cell_idx = findall(cell_arr[:,1] .!= 8)
        death_idx = 1
        for i = 1:n_cells
            # not taking jump instance into account right now
            # use so that I can increase it
            if cell_arr[i,1] == 0
                if 1-cdf(Gamma(1,1), time) < cell_arr[i,3]
                    # assign new cell state
                    cell_arr[i,1] = draw_fate()
                    # assign jump time
                    cell_arr[i,2] = time
                end
            end

            # check if cell state is not naive (0) or dead(8)
            # then check if cell will die
            if !(cell_arr[i,1] in (0,8))
                if (1-(cdf(Gamma(1,1), time-cell_arr[i,2]))) < cell_arr[i,4]
                    cell_arr[i,1] = 8
                    # assign jump time
                    cell_arr[i,2] = time
                end
            end

        end

        # get cell numbers
        n_th1[t] = sum(cell_arr[:,1] .== 1)
        n_tfh[t] = sum(cell_arr[:,1] .== 2)
        n_tr1[t] = sum(cell_arr[:,1] .== 3)
        n_naive[t] = sum(cell_arr[:,1] .== 0)    #println(t)
    end
    return [n_th1 n_tfh n_tr1]
end
# check how update rules for cells will apply
# find old python script where I did this

res_arr = run_sim(10, n_cells, time_arr)


for i=1:10
    if i == 1
        display(plot(time_arr, res_arr[i], color = [:blue :red :grey]))
    else
        display(plot!(time_arr, res_arr[i], color = [:blue :red :grey]))
    end
end
