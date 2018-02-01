###############################################################################
#                                                                             #
#  This file is part of the julia module for Multi Objective Optimization     #
#  (c) Copyright 2017 by Aritra Pal, Hadi Charkhgard                          #
#                                                                             #
# This license is designed to guarantee freedom to share and change software  #
# for academic use, but restricting commercial firms from exploiting our      #
# knowhow for their benefit. The precise terms and conditions for using,      #
# copying, distribution, and modification follow. Permission is granted for   #
# academic research use. The license expires as soon as you are no longer a   # 
# member of an academic institution. For other uses, contact the authors for  #
# licensing options. Every publication and presentation for which work based  #
# on the Program or its output has been used must contain an appropriate      # 
# citation and acknowledgment of the authors of the Program.                  #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

#####################################################################
# Solution Poilishing using K-OPT                                   #
#####################################################################

@inbounds function generating_starting_solutions_for_k_opt(instance::Union{MOLPInstance, BOLPInstance}, model::CPLEX.Model, bin_var_ind::Vector{Int64}, starting_solutions::Union{Vector{MOPSolution}, Vector{BOPSolution}}, k_max::Int64, k_min::Int64, timelimit::Float64)
    t0 = time()
    if typeof(instance) == MOLPInstance
        new_infeasible_starting_solutions = MOPSolution[]
    else
        new_infeasible_starting_solutions = BOPSolution[]
    end
    lhs = zeros(1, length(starting_solutions[1].vars))
    starting_solutions = starting_solutions[shuffle([1:length(starting_solutions)...])]
    i = 1
    while i <= length(starting_solutions) && (time()-t0) <= timelimit
        lhs = zeros(1, length(starting_solutions[i].vars))
        for j in bin_var_ind
            if starting_solutions[i].vars[j] == 1.0
                lhs[1, j] = -1.0
            else
                lhs[1, j] = 1.0
            end
        end
        CPLEX.add_constrs!(model, lhs, '<', [float(k_max) - sum(starting_solutions[i].vars[bin_var_ind])])
        CPLEX.add_constrs!(model, lhs, '>', [float(k_min+1) - sum(starting_solutions[i].vars[bin_var_ind])])
        CPLEX.optimize!(model)
        try
            tmp = CPLEX.get_solution(model)
            if typeof(starting_solutions[1]) == MOPSolution
                tmp2 = MOPSolution(vars=tmp)
            end
            if typeof(starting_solutions[1]) == BOPSolution
                tmp2 = BOPSolution(vars=tmp)
            end
            compute_objective_function_value!(tmp2, instance)
            push!(new_infeasible_starting_solutions, tmp2)
        catch
        end
        del_constrs!(model, size(instance.A)[1]+1, size(instance.A)[1]+2)
        i += 1
    end 
    select_unique_sols(new_infeasible_starting_solutions[shuffle([1:length(new_infeasible_starting_solutions)...])])
end

@inbounds function K_OPT(instance::Union{MOLPInstance, BOLPInstance}, model::CPLEX.Model, bin_var_ind::Vector{Int64}, starting_solutions::Union{Vector{MOPSolution}, Vector{BOPSolution}}, k_max::Int64, k_min::Int64, params)
    t0 = time()
    if typeof(starting_solutions[1]) == MOPSolution
        tmp = MOPSolution[]
        p = size(instance.c)[1]
    end
    if typeof(starting_solutions[1]) == BOPSolution
        tmp = BOPSolution[]
        p = 2
    end
    sols_to_explore, non_dom_sols = starting_solutions, starting_solutions
    sols_to_explore = sols_to_explore[shuffle([1:length(sols_to_explore)...])]
    timelimit = params[:timelimit]
    while length(sols_to_explore) >= 1 && (time()-t0) <= params[:timelimit]
        params[:timelimit] = (timelimit - (time()-t0))/3
        t1::Float64 = time()
        if typeof(starting_solutions[1]) == MOPSolution
            tmp2 = instance.c[1, :]
            for k in 2:size(instance.c)[1]
                tmp2 += instance.c[k, :]
            end
        end
        if typeof(starting_solutions[1]) == BOPSolution
            tmp2 = instance.c1 + instance.c2
        end
        params[:timelimit] -= (time()-t0)
        CPLEX.set_obj!(model, tmp2)
        new_infeasible_starting_solutions = generating_starting_solutions_for_k_opt(instance, model, bin_var_ind, sols_to_explore, k_max, k_min, params[:timelimit])
        params[:timelimit] = timelimit - (time()-t0)
        if params[:timelimit] < (2*timelimit)/3
            params[:timelimit] = (2*timelimit)/3
        end
        if length(bin_var_ind) == size(instance.A)[2]
            tmp = FPH(instance, model, new_infeasible_starting_solutions, params)
        else
            tmp = FPH(instance, model, bin_var_ind, new_infeasible_starting_solutions, params)
        end
        if length(tmp) >= 1
            sols_to_explore, non_dom_sols = FPBH.return_queue_for_local_search_operators(tmp, non_dom_sols)
            sols_to_explore = sols_to_explore[shuffle([1:length(sols_to_explore)...])]
        else
            break
        end
    end
    params[:timelimit] = timelimit
    [non_dom_sols..., tmp...]
end

@inbounds function K_OPT(instance::Union{MOLPInstance, BOLPInstance}, bin_var_ind::Vector{Int64}, starting_solutions::Union{Vector{MOPSolution}, Vector{BOPSolution}}, k_max::Int64, k_min::Int64, params)
    model = cplex_model(instance, 1)
    K_OPT(instance, model, bin_var_ind, starting_solutions, k_max, k_min, params)
end

@inbounds function K_OPT(instance::Union{MOLPInstance, BOLPInstance}, model::CPLEX.Model, bin_var_ind::Vector{Int64}, starting_solutions::Union{Vector{MOPSolution}, Vector{BOPSolution}}, params)
    t0 = time()
    non_dom_sols = starting_solutions
    k = 2
    step = 1
    count = Int64[]
    timelimit = params[:timelimit]
    while k <= length(bin_var_ind) && (time()-t0) <= timelimit
        params[:timelimit] = timelimit - (time()-t0)
        non_dom_sols = K_OPT(instance, model, bin_var_ind, non_dom_sols, k, k-step, params)
        push!(count, length(non_dom_sols))
        if length(count) >= 2 && abs(count[end] - count[end-1]) <= k
            step += 1
        end
        k += (step + 1)
    end
    params[:timelimit] = timelimit
    non_dom_sols
end

@inbounds function K_OPT(instance::Union{MOLPInstance, BOLPInstance}, bin_var_ind::Vector{Int64}, starting_solutions::Union{Vector{MOPSolution}, Vector{BOPSolution}}, params)
    model = cplex_model(instance, 1)
    K_OPT(instance, model, bin_var_ind, starting_solutions, params)
end

#####################################################################
# Parallization of Solution Poilishing using K-OPT                  #
#####################################################################

@inbounds function Parallel_K_OPT(instance::Union{MOLPInstance, BOLPInstance}, bin_var_ind::Vector{Int64}, starting_solutions::Union{Vector{MOPSolution}, Vector{BOPSolution}}, params)
    procs_ = setdiff(procs(), myid())[1:params[:total_threads]]
    inds = Vector{Int64}[]
    for i in 1:params[:total_threads]
        if i < params[:total_threads]
        	push!(inds, [(i-1)*ceil(Int64, length(starting_solutions)/params[:total_threads])+1:(i*ceil(Int64, length(starting_solutions)/params[:total_threads]))...])
        else
        	push!(inds, [(i-1)*ceil(Int64, length(starting_solutions)/params[:total_threads])+1:length(starting_solutions)...])
      	end
  	end
    if typeof(instance) == MOLPInstance
    	non_dom_sols = Vector{Vector{MOPSolution}}(params[:total_threads])
    else
    	non_dom_sols = Vector{Vector{BOPSolution}}(params[:total_threads])
   	end
   	@sync begin
    	for i in 1:params[:total_threads]
      		@async begin
                non_dom_sols[i] = remotecall_fetch(K_OPT, procs_[i], instance, bin_var_ind, starting_solutions[inds[i]], params)
            end
        end
    end
    select_non_dom_sols(vcat(non_dom_sols...))
end

@inbounds function SOL_POL(instance::Union{MOLPInstance, BOLPInstance}, model::CPLEX.Model, bin_var_ind::Vector{Int64}, starting_solutions::Union{Vector{MOPSolution}, Vector{BOPSolution}}, params)
    if params[:total_threads] == 1 && !params[:parallelism]
        K_OPT(instance, model, bin_var_ind, starting_solutions, params)
    else
        Parallel_K_OPT(instance, bin_var_ind, sort_non_dom_sols(starting_solutions), params)
    end
end
