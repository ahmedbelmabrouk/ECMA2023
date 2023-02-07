include("data/22_ulysses_9.tsp")

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

function resolution_dualite(n,L,B,K,W_v,w_v,W,coordinates)
    l=distance(n,coordinates)
    m=Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_MIPDISPLAY" =>2,"CPX_PARAM_TILIM" => 900))
    @variable(m,alpha>=0)
    @variable(m,beta[i=1:n,j=1:n]>=0)
    @variable(m,x[i=1:n,j=1:n],Bin)
    @variable(m,y[i=1:n,k=1:K],Bin)
    @variable(m,gamma[k=1:K]>=0)
    @variable(m,phi[i=1:n,k=1:K]>=0)



    @constraint(m,[i=1:n,j=i:n,k=1:K;i!=j],x[i,j]+y[i,k]-y[j,k]<=1)
    @constraint(m,[i=1:n,j=1:n,k=1:K;i!=j],x[i,j]-y[i,k]+y[j,k]<=1)
    @constraint(m,[i=1:n,j=1:n,k=1:K;i!=j],-x[i,j]+y[i,k]+y[j,k]<=1)
    @constraint(m,[i=1:n,j=1:n;i!=j],alpha+beta[i,j]>=(lh[i]+lh[j])*x[i,j])
    @constraint(m,[k=1:K],sum(w_v[i]*y[i,k]+W_v[i]*phi[i,k] for i=1:n)+W*gamma[k]<=B)
    @constraint(m,[i=1:n,k=1:K],gamma[k]+phi[i,k]>=w_v[i]*y[i,k])
    @constraint(m,[i=1:n],sum(y[i,k] for k=1:K)==1)

    @objective(m,Min,sum(l[i,j]*x[i,j]+3*beta[i,j] for i=1:n,j=i+1:n)+L*alpha)
    optimize!(m)
    
    return sum(l[i,j]*value.(x)[i,j] for i =1:n, j=1:n),solve_time(m)
end
data=["data/10_ulysses_3.tsp","data/10_ulysses_6.tsp","data/10_ulysses_9.tsp","data/14_burma_3.tsp","data/14_burma_6.tsp","data/14_burma_9.tsp","data/22_ulysses_3.tsp","data/22_ulysses_6.tsp","data/22_ulysses_9.tsp","data/26_eil_3.tsp"]



for path in data 
    open("results_dualité.txt","a") do file 
        include(path)
        ms=resolution_dualite(n,L,B,K,W_v,w_v,W,coordinates)
        println(file,string(n)*"&"*string(K)*"&"*string(ms[1])*"&"*string(ms[2]))
    end
end