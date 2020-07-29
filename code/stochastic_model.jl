##################################### add if I want module module stochastic_model

using Distributions
using HDF5
# param dict


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

# cell index dict
function draw_fate(d_fate, p1, p2, p3)
    # draw random number and compare to probabilities then choose fate
    n = rand(Float64)
    if 0 <= n < p1
        return d_fate["Th1"]
    elseif p1 <= n < p1+p2
        return d_fate["Tfh"]
    else
        return d_fate["Tr1"]
    end
end


function pos_fb(c, EC50)
    return c/(c+EC50)
end


function get_myc(time, deg_myc)
    return exp(-deg_myc*time)
end
# for each time point loop over each cell

function naive_diff(cell_arr, i, time, d_param, d_state, d_fate)
# note that naive cell does not need to take jump instance into account because last_state change = 0
    if 1-cdf(Gamma(d_param["alpha"],1/d_param["beta"]), time) < cell_arr[i, d_state["prob_change"]]
        # assign new cell state
        cell_arr[i,d_state["state"]] = d_fate["Prec"]

        # assign jump time
        # add this only because only effecotr cells die
        cell_arr[i, [d_state["last_change"], d_state["last_div"], d_state["last_death"]]] .= time
        cell_arr[i, d_state["prob_change"]] = rand(Float64)
    end
end


function prec_diff(cell_arr, i, time, d_param, d_state, d_fate, shift_arr,
    prob_arr)
# note that naive cell does not need to take jump instance into account because last_state change = 0
    if 1-cdf(Gamma(d_param["alpha"],1/d_param["beta"]), time-cell_arr[i, d_state["last_change"]]) <
        cell_arr[i, d_state["prob_change"]]
        # assign new cell state
        # unpack prob arr probabilities in draw fate
        fate = draw_fate(d_fate, prob_arr[i, :]...)
        cell_arr[i, d_state["state"]] = fate
        # update fate for gene shift
        shift_arr[i] = fate
        # assign jump time
        cell_arr[i, [d_state["last_change"], d_state["last_div"], d_state["last_death"]]] .= time
        cell_arr[i, d_state["prob_change"]] = rand(Float64)
    end
end


function cell_prolif(cell_arr, i, time, t_new, d_param, d_state, d_fate, mother_arr, daughter_arr)

    if cell_arr[i, d_state["prob_div0"]] > cell_arr[i, d_state["prob_div"]]

        create_cell(cell_arr, i, time, d_state, d_fate, mother_arr, daughter_arr)

        # for mother cell update division time and assign new division probability
        cell_arr[i, d_state["prob_div"]] = rand(Float64)
        cell_arr[i, d_state["prob_div0"]] = 0
        cell_arr[i, d_state["last_div"]] = time
    else
        cell_arr[i, d_state["prob_div0"]] = cell_arr[i, d_state["prob_div0"]] +
        cdf(Exponential(1/d_param["r_div"]), t_new-cell_arr[i,d_state["last_div"]]) -
        cdf(Exponential(1/d_param["r_div"]), time-cell_arr[i,d_state["last_div"]])
    end
end

function create_cell(cell_arr, i, time, d_state, d_fate, mother_arr, daughter_arr)
    # find a dead cell and make it alive
    # from all dead cells start with the first index then increase
    k = findfirst(cell_arr[:, d_state["state"]] .== d_fate["Dead"])
    cell_arr[k, d_state["state"]] = cell_arr[i, d_state["state"]]
    # for daughter cell add new div and death probabilities
    # also assign new last death/div/state change time
    cell_arr[k, d_state["prob_div"]] = rand(Float64)
    cell_arr[k, d_state["prob_death"]] = rand(Float64)
    # no one cares about state change at this point but anyway
    cell_arr[k, d_state["last_change"]] = time
    cell_arr[k, d_state["last_div"]] = time
    cell_arr[k, d_state["last_death"]] = time
    cell_arr[k, d_state["prob_div0"]] = 0

    # indicate that daugther cell genes will be subject to noise
    daughter_arr[k] = true
    mother_arr[i] = true
end


function cell_death(cell_arr, i, time, d_param, d_state, d_fate)

    if (1-(cdf(Exponential(1/d_param["r_death"]), time-cell_arr[i,d_state["last_death"]]))) <
        cell_arr[i,d_state["prob_death"]]
        cell_arr[i,d_state["state"]] = d_fate["Dead"]
    end
end


function update_cells(cell_arr, time, t_new, d_param, d_state,
    d_fate, shift_arr, mother_arr, daughter_arr, prob_arr)

    alive_cells = findall(cell_arr[:, d_state["state"]] .!= d_fate["Dead"])

    for i in alive_cells
        # note that the order of death vs prolif matters here
        # if they cells die first they cannot proliferate, reverse is not true
        # still not sure if I should use update state at t+1 so that all effects
        # differentiation
        if cell_arr[i, d_state["state"]] == d_fate["Naive"]
            naive_diff(cell_arr, i, time, d_param, d_state, d_fate)
        elseif cell_arr[i, d_state["state"]] == d_fate["Prec"]
            prec_diff(cell_arr, i, time, d_param, d_state, d_fate, shift_arr,
            prob_arr)
        else
            cell_prolif(cell_arr, i, time, t_new, d_param, d_state, d_fate,
            mother_arr, daughter_arr)
            cell_death(cell_arr, i, time, d_param, d_state, d_fate)
        end
    end
end


function update_cell_numbers(cell_arr, t, n_th1, n_tfh, n_tr1, n_naive, n_prec,
    d_state, d_fate)
    n_th1[t] = sum(cell_arr[:, d_state["state"]] .== d_fate["Th1"])
    n_tfh[t] = sum(cell_arr[:,d_state["state"]] .== d_fate["Tfh"])
    n_tr1[t] = sum(cell_arr[:, d_state["state"]] .== d_fate["Tr1"])
    n_naive[t] = sum(cell_arr[:, d_state["state"]] .== d_fate["Naive"])    #println(t)
    n_prec[t] = sum(cell_arr[:, d_state["state"]] .== d_fate["Prec"])
end


function create_cell_arr(n_cells, n_dead_cells, n_states, d_state, d_fate)
    cell_arr = zeros((n_cells, n_states))
    # assign random numbers for naive prec. transition
    cell_arr[:, [d_state["prob_change"], d_state["prob_death"], d_state["prob_div"]]] = rand(n_cells, 3)
    # assign random numbers for death event of cells
    dead_cell_arr = zeros((n_dead_cells, n_states))
    dead_cell_arr[:, d_state["state"]] .= d_fate["Dead"]
    cell_arr = [cell_arr; dead_cell_arr]

    return cell_arr
end

function gene_noise(gene_arr, mother_arr, daughter_arr, noise = 0.05)
    # for all mother cells take expr vals for all genes and assign to daughter cells with noise
    gene_arr[daughter_arr, :] = [rand(Normal(x, noise)) for x in gene_arr[mother_arr,:]]
end


function gene_shift(gene_arr, cells_idc, fate_ids)
    # from all cells that change fate take those genes that belong to new fate
    # update the gene expression (new val = old val + (1-old val) / 2)
    gene_arr[cells_idc, fate_ids] = gene_arr[cells_idc, fate_ids] .+
    ((1 .- gene_arr[cells_idc, fate_ids]) ./ 2)
end


function update_genes(gene_arr, shift_arr, mother_arr, daughter_arr)
    # find th1, tfh, tr1 idx in shift_arr and subset alive  cells
    # note that shift arr indices were updated only with respect to alive_cells
    th1_shift = shift_arr .== 2
    tfh_shift = shift_arr .== 3
    tr1_shift = shift_arr .== 4
    fate_th1 = 1:20
    fate_tfh = 21:40
    fate_tr1 = 41:60
    # update gene shift for th1, tfh and tr1 cells
    gene_shift(gene_arr, th1_shift, fate_th1)
    gene_shift(gene_arr, tfh_shift, fate_tfh)
    gene_shift(gene_arr, tr1_shift, fate_tr1)
    # add noise to daughter cells based on prolif arr
    gene_noise(gene_arr, mother_arr, daughter_arr)
end


function create_prob_arr(gene_arr, fate_th1, fate_tfh, fate_tr1, d_param, mode = "fixed")
    # take gene arr and collapse into probability array for fate decision
    # fate_th1... should be range objects
    # returns prob arr, same n row as gene_arr with 3 cols

    if mode == "predetermined"
        fate_idx = [fate_th1, fate_tfh, fate_tr1]
        # take the sum across each set of signature genes
        prob_arr = [sum(gene_arr[:,idx], dims=2) for idx in fate_idx]
        # concatenate into one horizontal vector
        prob_arr = hcat(prob_arr...)
        # normalize the probability
        prob_arr = prob_arr ./ sum(prob_arr, dims=2)
    else
        # create prob arr based on param dict
        # make arr of ones and multiply by column with default probability
        n, m = size(gene_arr)
        prob_arr = ones(n, 3)
        probs = ["prob_Th1", "prob_Tfh", "prob_Tr1"]
        for (index, value) in enumerate(probs)
            prob_arr[:,index] = prob_arr[:,index] .* d_param[value]
        end

    end

    return prob_arr
end

function stoc_model(n_cells, n_genes, time_arr, d_param, d_state, d_fate)
    ######################################################### create cell array
    n_dead_cells = 500
    n_states = length(d_state)
    cell_arr = create_cell_arr(n_cells, n_dead_cells, n_states, d_state, d_fate)
    gene_arr = rand(n_cells+n_dead_cells, n_genes)

    fate_th1 = 1:20
    fate_tfh = 21:40
    fate_tr1 = 41:60
    prob_arr = create_prob_arr(gene_arr, fate_th1, fate_tfh, fate_tr1, d_param)

    ################################################ create arr for cell numbers
    n_naive = zeros(Int64, length(time_arr))
    n_th1 = zeros(Int64, length(time_arr))
    n_tfh = zeros(Int64, length(time_arr))
    n_tr1 = zeros(Int64, length(time_arr))
    n_prec = zeros(Int64, length(time_arr))

    for t = 1:(length(time_arr)-1)
        time = time_arr[t]
        t_new = time_arr[t+1]

        #myc = get_myc(time, deg_myc)
        #r_div = r_div_base#*pos_fb(myc, EC50_myc)
        # find all alive cells, to only loop over them

        # make an array to remember which cells change fate for gene arr
        shift_arr = zeros(Int64, n_cells+n_dead_cells)
        # make an array to check which cells are newly created
        mother_arr = falses(n_cells+n_dead_cells)
        daughter_arr = falses(n_cells+n_dead_cells)

        update_cells(cell_arr, time, t_new, d_param, d_state,
        d_fate, shift_arr, mother_arr, daughter_arr, prob_arr)

        update_cell_numbers(cell_arr, t, n_th1, n_tfh, n_tr1, n_naive, n_prec,
        d_state, d_fate)

        update_genes(gene_arr, shift_arr, mother_arr, daughter_arr)
    end

    # set last element of cel_arr to timestep before because loop only goes so far
    cell_arr[end, :] = cell_arr[(end-1),:]
    update_cell_numbers(cell_arr, length(time_arr), n_th1, n_tfh, n_tr1, n_naive, n_prec,
    d_state, d_fate)
    # find alive cells
    alive_cells = findall(cell_arr[:, d_state["state"]] .!= d_fate["Dead"])
    # only choose genes of cells that are alive at the end of simulation
    gene_arr = gene_arr[alive_cells, :]
    cell_arr = cell_arr[alive_cells, d_state["state"]]
    # add ground truth fate to gene array
    gene_arr = [gene_arr cell_arr]

    df = [time_arr n_th1 n_tfh n_tr1 n_prec]
    return (df, gene_arr)
end

#################################################### change below if I want module
# second entry is jump time (last state transition)
# third entry is cumulative
function run_sim(n_sim, n_cells, n_genes, time_arr, d_param, d_state, d_fate)
    cell_arr = []
    gene_arr = []

    # run stoc_model several times, each time add gene arr and cell arr to list
    Threads.@threads for i = 1:n_sim
        res = stoc_model(n_cells, n_genes, time_arr, d_param, d_state, d_fate)
        cell_df = res[1]
        gene_df = res[2]
        push!(cell_arr, cell_df)
        push!(gene_arr, gene_df)
    end

    # concatenate cell arr to one large list
    res_cells = vcat(cell_arr...)

    return (res_cells, gene_arr)
end

# parameters
n_cells = 2000
# cell can be alive or dead
# create cell array
n_genes = 60
n_sim = 3

time_arr = range(0, 5, step = 0.001)
cell_arr, gene_arr = run_sim(n_sim, n_cells, n_genes, time_arr, d_param, d_state, d_fate)

# save files as hdf5 format
sc_data =h5open("home/burt/Documents/projects/2020/tcell_stochastic/output/scseq_sim.h5","w")
for i=1:n_sim
    sc_data[string(i)] = gene_arr[i]
end
close(sc_data)

cell_data =h5open("home/burt/Documents/projects/2020/tcell_stochastic/output/model_output.h5","w")
cell_data["cell_data"] = cell_arr
close(cell_data)
