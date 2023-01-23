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

function get_arc(sol::Matrix{Float64},n::Int64)
    path = Matrix{Int64}(zeros(n,2))
    cpt = 1
    for i in 1:n
        for j in 1:n
            if sol[i,j] == 1
                path[cpt,1] = i
                path[cpt,2] = j
                cpt = cpt + 1
                break
            end
        end
    end
    return path[1:(cpt-1),:]
end

function find_path(sol::Matrix{Float64},de::Int64,fi::Int64)
    current = de
    n = size(sol)[1]
    res = Vector{Int64}(zeros(120))
    res[1] = current
    cpt = 2
    for i in 1:n
        for j in 1:n
            if sol[current,j] == 1
                current = j
                res[cpt] = current
                cpt = cpt + 1
                break
            end
        end
        if current == fi
            break
        end
    end
    return res[1:(cpt-1)]
end

function is_int(v::Vector{Int64}, ind::Int64)
    for i in v
        if i == ind
            return 1
        end
    end
    return 0
end

function is_in_path(v1::Vector{Int64},v2::Vector{Int64})
    for i in v2
        for j in v1
            if i == j
                return 1
            end
        end
    end
    return 0
end

function DFJ_const(sol::Matrix{Float64},m::Matrix{Int64})
    res = Vector{Int64}(zeros(100))
    cpt = 1
    for a in m[:,1]
        if is_int(m[:,2],a) == 1
            path = find_path(sol,a,a)
            if path[1] == path[length(path)] && is_in_path(res,path) == 0
                res[cpt] = a
                cpt = cpt + 1
            end
        end
    end
    if cpt == 1
        return zeros(0)
    end
    if cpt > 1
        return res[1:(cpt-1)]
    end
end
