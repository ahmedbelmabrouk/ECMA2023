# PROJET ECMA : Ahmed Belmabrouk

using JuMP
using Dates
using CPLEX

include("data/10_ulysses_3.tsp")

function distance(n::Int64,coordinates::Matrix{Float64})
    l=Matrix{Float64}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            l[i,j]=sqrt((coordinates[i,1]-coordinates[j,1])^2+(coordinates[i,2]-coordinates[j,2])^2)
        end 
    end 
    return l
end
l=distance(n,coordinates)
donnee=n,L,B,K,W_v,w_v,W,coordinates,l

#**********Probleme maitre************
m = Model(CPLEX.Optimizer)

@variable(m, x[i in 1:n, j in 1:n] , Bin)
@variable(m,y[i=1:n,k=1:K],Bin)
@variable(m, z >=0)

#Définition des contraintes
@constraint(m, (sum(l[i,j] * x[i,j] for i in 1:n , j in i+1:n  )) <= z )
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K ;i!=j ],x[i,j] + y[i,k] - y[j,k] <= 1) #lien entre x et y 1/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j  ],x[i,j] - y[i,k] + y[j,k] <= 1) #lien entre x et y 2/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j  ],-x[i,j] + y[i,k] + y[j,k] <= 1) #lien entre x et y 3/3
@constraint(m,[k in 1:K], sum(w_v[i] * y[i,k] for i in 1:n ) <= B)
@constraint(m,[i in 1:n],sum(y[i,k] for k in 1:K ) == 1)
#@constraint(m,[i in 1:n],x[i,i] ==0)

#Définition de l'objectif
@objective(m, Min, z)
optimize!(m)
z_etoile=objective_value(m)
print("==> z* =",z_etoile,"<=======")
x_etoile=value.(x)
for i in 1:n 
    for j in 1:n 
        println("x[",i,",",j,"]",x_etoile[i,j])
    end
end
y_etoile=value.(y)