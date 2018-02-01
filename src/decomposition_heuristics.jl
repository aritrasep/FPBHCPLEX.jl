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
# Bi-Objective Pure and Mixed Binary Programs                       #
#####################################################################

#####################################################################
## Modified Perpendicular Search Method                            ##
#####################################################################

@inbounds function modified_perpendicular_search_method(instance::BOLPInstance, model::CPLEX.Model, bin_var_ind::Vector{Int64}, pts_to_explore::Vector{Vector{Float64}}, params)
    t0 = time()
    non_dom_sols = BOPSolution[]
    sols_to_explore = BOPSolution[]
    pos1, pos2 = length(instance.cons_lb), length(instance.cons_lb)
    timelimit = params[:timelimit]
    while length(pts_to_explore) >= 1 && time()-t0 <= timelimit
        params[:timelimit] = (timelimit - (time()-t0))/length(pts_to_explore)
        if length(pts_to_explore) == 1 && pts_to_explore[1] == [Inf, Inf]
            params[:timelimit] = params[:timelimit] / 10.0
        end
        current_pt_to_explore = splice!(pts_to_explore, 1)
        for i in 1:2
            if current_pt_to_explore[i] != Inf
                if i == 1
                    CPLEX.add_constr!(model, -1.0*instance.c1, '>', -1.0*current_pt_to_explore[1])
                    instance.A = vcat(instance.A, -1.0*instance.c1')
                else
                    CPLEX.add_constr!(model, -1.0*instance.c2, '>', -1.0*current_pt_to_explore[2])
                    instance.A = vcat(instance.A, -1.0*instance.c2')
                end
                instance.cons_lb = vcat(instance.cons_lb, -1.0*current_pt_to_explore[i])
                instance.cons_ub = vcat(instance.cons_ub, Inf)
                pos2 += 1
            end
        end
        tmp = FPH(instance, model, bin_var_ind, params)
        if length(tmp) >= 1
            push!(sols_to_explore, tmp...)
        end
        if pos1 < pos2
            instance.A = instance.A[1:pos1,:]
            instance.cons_lb = instance.cons_lb[1:pos1]
            instance.cons_ub = instance.cons_ub[1:pos1]
            del_constrs!(model, pos1+1, pos2)
            pos2 = pos1
        end
        if length(pts_to_explore) == 0 && length(sols_to_explore) >= 1
            pts_to_explore, non_dom_sols = FPBH.return_queue_for_modified_perpendicular_search_method(sols_to_explore, non_dom_sols, (timelimit - (time()-t0))/3)
            sols_to_explore = BOPSolution[]
        end
    end
    params[:timelimit] = timelimit
    [non_dom_sols..., sols_to_explore...]
end

#####################################################################
## Parallel Modified Perpendicular Search Method                   ##
#####################################################################

@inbounds function parallel_modified_perpendicular_search_method(instance::BOLPInstance, bin_var_ind::Vector{Int64}, params)
    t0 = time()
    pts_to_explore = [[Inf, Inf]]
    non_dom_sols = BOPSolution[]
    sols_to_explore = BOPSolution[]
    procs_ = setdiff(procs(), myid())[1:params[:total_threads]]
    timelimit = params[:timelimit]
    tmp = fill(BOPSolution[], length(procs_))
    while length(pts_to_explore) >= 1 && time()-t0 <= timelimit
        params[:timelimit] = (timelimit - (time()-t0))/length(pts_to_explore)
        params[:timelimit] = params[:timelimit] * params[:total_threads]    
        if length(pts_to_explore) == 1 && pts_to_explore[1] == [Inf, Inf]
            params[:timelimit] = params[:timelimit] / 10.0
        end
        @sync begin
            for i in 1:length(procs_)
                @async begin
                    if length(pts_to_explore) >= 1
                        tmp[i] = remotecall_fetch(FPH, procs_[i], copy(instance), bin_var_ind, splice!(pts_to_explore, 1), params)
                    end
                end
            end
        end
        for i in 1:length(procs_)
            if length(tmp[i]) >= 1
                push!(sols_to_explore, tmp[i]...)
            end
        end
        if length(pts_to_explore) == 0 && length(sols_to_explore) >= 1
            pts_to_explore, non_dom_sols = FPBH.return_queue_for_modified_perpendicular_search_method(sols_to_explore, non_dom_sols, (timelimit - (time()-t0))/3)
            sols_to_explore = BOPSolution[]
        end
    end
    [non_dom_sols..., sols_to_explore...]
end

@inbounds function modified_perpendicular_search_method(instance::BOLPInstance, model::CPLEX.Model, bin_var_ind::Vector{Int64}, params)
    pts_to_explore = [[Inf, Inf]]
    modified_perpendicular_search_method(instance, model, bin_var_ind, pts_to_explore, params)
end

@inbounds function MPSM(instance::BOLPInstance, model::CPLEX.Model, bin_var_ind::Vector{Int64}, params)
    if params[:total_threads] == 1 && !params[:parallelism]
        modified_perpendicular_search_method(instance, model, bin_var_ind, params)
    else
        parallel_modified_perpendicular_search_method(instance, bin_var_ind, params)
    end
end

#####################################################################
# Multiobjective Pure or Mixed Binary Programs                      #
#####################################################################

#####################################################################
## Modified Full P Split Method                                    ##
#####################################################################

@inbounds function modified_full_p_split_method(instance::MOLPInstance, model::CPLEX.Model, bin_var_ind::Vector{Int64}, pts_to_explore::Vector{Vector{Float64}}, params)
    t0 = time()
    p = size(instance.c)[1]
    non_dom_sols = MOPSolution[]
    sols_to_explore = MOPSolution[]
    pos1, pos2 = length(instance.cons_lb), length(instance.cons_lb)
    timelimit = params[:timelimit]
    while length(pts_to_explore) >= 1 && time()-t0 <= timelimit
        params[:timelimit] = (timelimit - (time()-t0))/length(pts_to_explore)
        if length(pts_to_explore) == 1 && pts_to_explore[1] == fill(Inf, size(instance.c)[1])
            params[:timelimit] = params[:timelimit] / 10.0
        end
        current_pt_to_explore = splice!(pts_to_explore, 1)
        for i in 1:p
            if current_pt_to_explore[i] != Inf
                CPLEX.add_constr!(model, -1.0*vec(instance.c[i, :]), '>', -1.0*current_pt_to_explore[i])
                instance.A = vcat(instance.A, -1.0*instance.c[i, :]')
                instance.cons_lb = vcat(instance.cons_lb, -1.0*current_pt_to_explore[i])
                instance.cons_ub = vcat(instance.cons_ub, Inf)
                pos2 += 1
            end
        end
        tmp = FPH(instance, model, bin_var_ind, params)
        if length(tmp) >= 1
            push!(sols_to_explore, tmp...)
        end
        if pos1 < pos2
            instance.A = instance.A[1:pos1,:]
            instance.cons_lb = instance.cons_lb[1:pos1]
            instance.cons_ub = instance.cons_ub[1:pos1]
            del_constrs!(model, pos1+1, pos2)
            pos2 = pos1
        end
        if length(pts_to_explore) == 0 && length(sols_to_explore) >= 1
            pts_to_explore, non_dom_sols = FPBH.return_queue_for_modified_full_p_split_method(sols_to_explore, non_dom_sols, (timelimit - (time()-t0))/3)
            sols_to_explore = MOPSolution[]
        end
    end
    params[:timelimit] = timelimit
    if length(non_dom_sols) >= 1 || length(sols_to_explore) >= 1
        select_non_dom_sols([non_dom_sols..., sols_to_explore...])
    else
        non_dom_sols
    end
end

#####################################################################
## Parallel Modified P Split Search Method                         ##
#####################################################################

@inbounds function parallel_modified_full_p_split_method(instance::MOLPInstance, bin_var_ind::Vector{Int64}, params)
    t0 = time()
    pts_to_explore = [fill(Inf, size(instance.c)[1])]
    non_dom_sols = MOPSolution[]
    sols_to_explore = MOPSolution[]
    procs_ = setdiff(procs(), myid())[1:params[:total_threads]]
    timelimit = params[:timelimit]
    tmp = fill(MOPSolution[], length(procs_))
    while length(pts_to_explore) >= 1 && time()-t0 <= timelimit
        params[:timelimit] = (timelimit - (time()-t0))/length(pts_to_explore)
        params[:timelimit] = params[:timelimit] * params[:total_threads]    
        if length(pts_to_explore) == 1 && pts_to_explore[1] == fill(Inf, size(instance.c)[1])
            params[:timelimit] = params[:timelimit] / 10.0
        end
        @sync begin
            for i in 1:length(procs_)
                @async begin
                    if length(pts_to_explore) >= 1
                        tmp[i] = remotecall_fetch(FPH, procs_[i], copy(instance), bin_var_ind, splice!(pts_to_explore, 1), params)
                    end
                end
            end
        end
        for i in 1:length(procs_)
            if length(tmp[i]) >= 1
                push!(sols_to_explore, tmp[i]...)
            end
        end
        if length(pts_to_explore) == 0 && length(sols_to_explore) >= 1
            pts_to_explore, non_dom_sols = FPBH.return_queue_for_modified_full_p_split_method(sols_to_explore, non_dom_sols, (timelimit - (time()-t0))/3)
            sols_to_explore = MOPSolution[]
        end
    end
    params[:timelimit] = timelimit
    if length(non_dom_sols) >= 1 || length(sols_to_explore) >= 1
        select_non_dom_sols([non_dom_sols..., sols_to_explore...])
    else
        non_dom_sols
    end
end

@inbounds function modified_full_p_split_method(instance::MOLPInstance, model::CPLEX.Model, bin_var_ind::Vector{Int64}, params)
    pts_to_explore::Vector{Vector{Float64}} = [fill(Inf, size(instance.c)[1])]
    modified_full_p_split_method(instance, model, bin_var_ind, pts_to_explore, params)
end

@inbounds function MFPSM(instance::MOLPInstance, model::CPLEX.Model, bin_var_ind::Vector{Int64}, params)
    if params[:total_threads] == 1 && !params[:parallelism]
        modified_full_p_split_method(instance, model, bin_var_ind, params)
    else
        parallel_modified_full_p_split_method(instance, bin_var_ind, params)
    end
end
