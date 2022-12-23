using JuMP, Cbc, DelimitedFiles, CSV, DataFrames
function f()

    dat = readdlm("donor-pool2.csv", ',', '\n', comments=true)

    fr = dat[:,1]
    to = dat[:,2]
    w = dat[:,3]
    N = dat[:,4]
    N = N[1:3,]

    V = union(fr,to)       # set of all nodes
    E = collect(zip(fr,to)) # set of all edges

    W = Dict( (i,j) => k for (i,j,k) in zip(fr,to,w) )  # weights on edges

    m = Model()
    set_optimizer(m, Cbc.Optimizer)

    @variable(m, y[(i,j) in E], Bin)

    @variable(m, f_in[i in V])

    @variable(m, f_out[i in V])

    @objective(m, Max, sum( W[(i,j)]*y[(i,j)] for (i,j) in E ))

    # 1
    @constraint(m, flow_out[k in V], sum(y[(i,j)] for (i,j) in E if i==k) == f_out[k])
    @constraint(m, flow_in[k in V], sum(y[(i,j)] for (i,j) in E if j==k) == f_in[k])

    # 2
    @constraint(m, flow_ineq[k in setdiff(V,N)], f_out[k] <= f_in[k])
    @constraint(m, flow_ineq2[k in setdiff(V,N)], f_in[k] <= 1)

    # 3
    @constraint(m, flow_out_ndd[k in N], f_out[k] <= 1)

    Cycle = []
    # Note that keeping the list of subtours in Subtour is not really
    # necessary. I just want to see the subtours that are identified so
    # that I know the algorithm is working.

    # loop to iteratively identify subtours and add the violated
    # constraints to the model.
    while true # all cycles aren't under 5
        cycles_over_five = 0
        JuMP.optimize!(m)
        set = []
        for (i,j) in E
            if value.(y[(i,j)]) > 0
                set = union(set,i,j)
            end
        end
        print(length(set))
        #print(value.(y))
        # get the pairs that are a part of the solution
        searched = []
        U = []
        ExtendU = set[rand(1:end)] # choose a node at random
        push!(searched,ExtendU)
        #print("U ",ExtendU,"\n")


        # look for the next edge in our solution to
        # map all pairs to one another. capture whether
        # the pairs are part of a cycle in which we will
        # record that cycle, or if they are a part of a
        # chain in which we will ignore.
        while length(setdiff(set,searched))!=0 # while we haven't searched all edges
            push!(U,ExtendU)
            ExtendU = nothing
            for j in union(U[1],setdiff(set,searched))
                if (in((U[end],j),E) && JuMP.value(y[(U[end],j)]) > 0)
                    ExtendU = j
                    print(length(searched),"\n")
                    push!(searched,ExtendU)
                    break
                end
            end
            # if we've found the end of a cycle
            if U[1]==ExtendU
                print("Cycle Found","\n")
                # we've found a cycle, lets add it to our cycle list and create
                # a constraint if necessary
                #print("hi")
                if length(U) > 5
                    push!(Cycle,U)
                    cycles_over_five = cycles_over_five + 1
                end
                U = []
                rn = setdiff(set,searched)
                ExtendU = rn[rand(1:end)]
                push!(searched,ExtendU)
            end
            # if we've found the end of a chain
            if isnothing(ExtendU)
                print("Chain Found","\n")
                U = []
                rn = setdiff(set,searched)
                ExtendU = rn[rand(1:end)] # don't pick the same one twice..
                push!(searched,ExtendU)
                # we've found a chain!
            end
        end
        for z in 1:length(Cycle)
            print("constraint added")
            @constraint(m, sum(y[(i,j)] for i in Cycle[z], j in Cycle[z] if in((i,j),E)) <= length(Cycle[z]) - 1)
        end
        if cycles_over_five==0
            JuMP.optimize!(m)
            for (i,j) in E
                if value.(y[(i,j)]) > 0
                    println(i, " to ", j)
                end
            end
            break
        end
    end
    return(JuMP.value.(y),E)
end

function find_chains(y,E,soln,ndd)
    Chain = []
    for i in ndd
        stop=false
        searched = []
        U = []
        ExtendU = i # start at ndd to find chain
        push!(searched,ExtendU)

        while stop==false
        # convert solution to a usable format for health professionals
            push!(U,ExtendU)
            ExtendU = nothing
            for j in union(U[1],setdiff(soln,searched))
                if (in((U[end],j),E) && y[(U[end],j)] > 0)
                    ExtendU = j
                    print(length(searched),"\n")
                    push!(searched,ExtendU)
                    break
                end
            end
            if isnothing(ExtendU)
                print("Chain Found","\n")
                push!(Chain,U)
                U = []
                rn = setdiff(soln,searched)
                ExtendU = rn[rand(1:end)] # don't pick the same one twice..
                push!(searched,ExtendU)
                stop=true
                # we've found a chain!
            end
        end
    end
    return(Chain)
end
# find cycles
function find_cycles(foundchains,y,E,soln)
    fcarray = []
    Cycle = []
    for i in 1:length(foundchains)
        union!(fcarray,foundchains[i])
    end
    tosearch = setdiff(soln,fcarray) # search only pairs a part of cycles
    searched = []
    U = []
    ExtendU = tosearch[rand(1:end)] # start at ndd to find chain
    push!(searched,ExtendU)
    while length(setdiff(tosearch,searched))!=0 # while we haven't searched all edges
        push!(U,ExtendU)
        ExtendU = nothing
        for j in tosearch
            if (in((U[end],j),E) && y[(U[end],j)] > 0)
                ExtendU = j
                print(length(searched),"\n")
                push!(searched,ExtendU)
                break
            end
        end
        if U[1]==ExtendU
            # found a cycle
            push!(U,ExtendU)
            push!(Cycle,U)
            U = []
            rn = setdiff(tosearch,searched)
            ExtendU = rn[rand(1:end)]
            push!(searched,ExtendU)
        end
    end
    push!(U,ExtendU)
    push!(Cycle,U)
    return(Cycle)
end
y, E = f()
soln = []
for (i,j) in E
    if y[(i,j)] > 0
        union!(soln,i,j)
    end
end
ndd = [13,978,999] # enter the NDDs here
foundchains = find_chains(y,E,soln,ndd)
Cycle = find_cycles(foundchains,y,E,soln)

CSV.write("cycles.csv",DataFrame(Cycle),delim=',')
CSV.write("chains.csv",DataFrame(foundchains),delim=',')
