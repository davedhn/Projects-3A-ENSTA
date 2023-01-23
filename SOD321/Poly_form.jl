using JuMP, Gurobi

function dif(index::Int64, len::Int64)
    res = Vector{Int64}(zeros(n-1))
    c = 1
    for i in 1:n
        if i!=index
            res[c] = i
            c = c + 1
        end
    end
    return res
end

function ind_flow(i::Int64)
    if i == d
        return 1
    end
    if i == f
        return -1
    end
    if i != d && i != f
        return 0
    end
end

function Prob_w_poly_const(n::Int64,d::Int64,f::Int64,Amin::Int64,Nr::Int64,R::Int64,regions::Dict,D::Matrix{Int64})
    # Model
    model = Model(Gurobi.Optimizer)

    # Variables
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, u[1:n], Int)

    # Objective
    @objective(model, Min, sum(sum(D[i,j]*x[i,j] for j in dif(i,n)) for i in 1:n))

    # Starts at d
    @constraint(model, sum(x[d,j] for j in dif(d,n)) == 1)
    @constraint(model, sum(x[i,d] for i in dif(d,n)) == 0)

    # Ends at f
    @constraint(model, sum(x[i,f] for i in dif(f,n)) == 1)
    @constraint(model, sum(x[f,j] for j in dif(f,n)) == 0)

    # Come from no more than one city
    for j in 1:n
        @constraint(model, sum(x[i,j] for i in dif(j,n)) <= 1)
    end

    for i in 1:n
        @constraint(model, sum(x[i,j] for j in dif(i,n)) <= 1)
    end

    # Flow conservation
    for i in 1:n
        @constraint(model, sum(x[i,j] for j in dif(i,n)) - sum(x[j,i] for j in dif(i,n)) == ind_flow(i))
    end

    # Min number of city to visit
    @constraint(model, sum(sum(x[i,j] for j in dif(i,n)) for i in 1:n) >= Amin - 1)

    # Max distance to land
    for i in 1:n
        @constraint(model, sum(x[i,j]*D[i,j] for j in dif(i,n)) <= R)
    end

    # All regions must be visited
    for num_reg in 1:Nr
        @constraint(model, sum(sum(x[i,j]+x[j,i] for j in dif(i,n)) for i in regions[num_reg]) >= 1)
    end

    #=
    # Only a single tour covering all cities (MTZ constraints)
    for i in 1:n
        for j in 1:n
            if i!=d && j!=f
                @constraint(model, u[i]+1-n*(1-x[i,j])<=u[j])
            end
        end
    end

    for i in 1:n
        @constraint(model, u[i]>=1)
        @constraint(model, u[i]<=n-1)
    end

    @constraint(model, u[d] == 1)
    @constraint(model, u[f] == n-1)
    =#

    #=
    # Only a single tour covering all cities (SF constraints)
    @variable(model, q[1:n,1:n])
    @variable(model, z[1:n], Bin)

    for i in 1:n
        for j in 1:n
            @constraint(model, q[i,j] >= 0)
        end
    end


    for i in 1:n
        for j in 1:n
            @constraint(model, q[i,j]<=(n-1)*x[i,j])
        end
    end

    @constraint(model, sum(q[d,j] for j in 1:n) == sum(z[k] for k in dif(d,n)))

    for k in dif(d,n)
        @constraint(model, sum(q[i,k] for i in 1:n) - sum(q[k,j] for j in 1:n) == z[k])
    end

    for k in dif(d,n)
        @constraint(model, sum(x[i,k] for i in 1:n) == z[k])
    end
    =#


    # Only a single tour covering all cities (MCF constraints)
    @variable(model, q[1:n,1:n,1:n])
    @variable(model, z[1:n], Bin)

    for i in 1:n
        for j in 1:n
            for k in dif(d,n)
                @constraint(model, q[i,j,k] >= 0)
            end
        end
    end


    for i in 1:n
        for j in 1:n
            for k in dif(d,n)
                @constraint(model, q[i,j,k] <=x[i,j])
            end
        end
    end

    @constraint(model, sum(x[d,j] for j in 1:n) == 1)
    @constraint(model, sum(x[j,f] for j in 1:n) == 1)

    for i in 1:n
        for k in dif(d,n)
            if i == d
                @constraint(model, sum(q[i,j,k] for j in 1:n) - sum(q[j,i,k] for j in 1:n) == z[k])
            end
            if i == k
                @constraint(model, sum(q[i,j,k] for j in 1:n) - sum(q[j,i,k] for j in 1:n) == -z[k])
            end
            if i != d && i != k
                @constraint(model, sum(q[i,j,k] for j in 1:n) - sum(q[j,i,k] for j in 1:n) == 0)
            end
        end
    end

    for k in dif(d,n)
        @constraint(model, sum(x[i,k] for i in 1:n) == z[k])
    end
    

    # Resolution
    JuMP.optimize!(model)

    # Print results
    obj_val = JuMP.objective_value(model)
    println("Objective value : ", obj_val)

    sol = JuMP.value.(x)

    path = string(d)
    current = d
    for i in 1:n
        if current == f
            break
        end
        for j in dif(d,n)
            if sol[current,j] == 1
                path = path*"->"*string(j)
                current = j
                break
            end
        end
    end
    println("The path is :", path)
end
