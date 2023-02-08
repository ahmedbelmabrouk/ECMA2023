

using JuMP
using CPLEX
function distance(n::Int64,coordinates::Matrix{Float64})
    l=Matrix{Float64}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            l[i,j]=sqrt((coordinates[i,1]-coordinates[j,1])^2+(coordinates[i,2]-coordinates[j,2])^2)
        end 
    end 
    return l
end

function resolution_statique(n,L,B,K,W_v,w_v,W,coordinates)
    
    l=distance(n,coordinates)
    m=Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_MIPDISPLAY" =>2,"CPX_PARAM_TILIM" => 900))
    @variable(m,x[i=1:n,j=1:n],Bin)
    @variable(m,y[i=1:n,k=1:K],Bin)
    @variable(m, z >=0)
    @constraint(m, (sum(l[i,j] * x[i,j] for i in 1:n , j in i+1:n  )) <= z )
    @constraint(m,[i in 1:n ,j in 1:n , k in 1:K ;i!=j ],x[i,j] + y[i,k] - y[j,k] <= 1) #lien entre x et y 1/3
    @constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j  ],x[i,j] - y[i,k] + y[j,k] <= 1) #lien entre x et y 2/3
    @constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j  ],-x[i,j] + y[i,k] + y[j,k] <= 1) #lien entre x et y 3/3
    @constraint(m,[k in 1:K], sum(w_v[i] * y[i,k] for i in 1:n ) <= B)
    @constraint(m,[i in 1:n],sum(y[i,k] for k in 1:K ) == 1)
    @objective(m,Min,z)
    return m 
end
data=["data/10_ulysses_3.tsp","data/10_ulysses_6.tsp","data/10_ulysses_9.tsp","data/14_burma_3.tsp","data/14_burma_6.tsp","data/14_burma_9.tsp","data/22_ulysses_3.tsp","data/22_ulysses_6.tsp","data/22_ulysses_9.tsp","data/26_eil_3.tsp"]
for path in data 
    open("results_statique.txt","a") do file 
        include(path)
        ms=resolution_statique(n,L,B,K,W_v,w_v,W,coordinates)
        optimize!(ms)
        println(file,string(n)*"&"*string(K)*"&"*string(objective_value(ms))*"&"*string(solve_time(ms)))
    end
end