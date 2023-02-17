# PROJET ECMA : Ahmed Belmabrouk & Rym Ben Maktouf

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

#***********Max des poids (poids au pire des cas )************
w_v_val=[w_v[i] *(1+ W_v[i]) for i in 1:n]

#***************Heuristique****************************
#**********Probleme maitre************
m = Model(CPLEX.Optimizer)

@variable(m, x[i in 1:n, j in 1:n] , Bin)
@variable(m,y[i=1:n,k=1:K],Bin)
@variable(m, z >=0)

#Définition des contraintes
@constraint(m, (sum(l[i,j] * x[i,j] for i in 1:n , j in i+1:n )) <= z )
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j ],x[i,j] + y[i,k] - y[j,k] <= 1) #lien entre x et y 1/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j ],x[i,j] - y[i,k] + y[j,k] <= 1) #lien entre x et y 2/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j ],-x[i,j] + y[i,k] + y[j,k] <= 1) #lien entre x et y 3/3
@constraint(m,[k in 1:K], sum(w_v_val[i] * y[i,k] for i in 1:n ) <= B)
@constraint(m,[i in 1:n],sum(y[i,k] for k in 1:K ) == 1)

#Définition de l'objectif
@objective(m, Min, z)
optimize!(m)
z_etoile=objective_value(m)
print("==> z* =",z_etoile,"<=======")
x_etoile=value.(x)
y_etoile=value.(y)
#**********sous probleme 1************
m1 = Model(CPLEX.Optimizer)

@variable(m1, sigma1[i in 1:n,j in 1:n] >=0)
@constraint(m1, sum(sigma1[i,j] for i in 1:n,j in i+1:n) <= L)
@constraint(m1 , [i in 1:n, j in 1:n] , sigma1[i,j] <= 3)
@objective(m1, Max, sum((l[i,j]+sigma1[i,j]*(lh[i]+lh[j])) * x_etoile[i,j] for i in 1:n, j in i+1:n))
optimize!(m1)
z1=objective_value(m1)
vsigma1=value.(sigma1)
while (z1 - z_etoile) > 1e-4 
    @constraint(m, (sum((l[i,j]+ vsigma1[i,j]*(lh[i]+lh[j]))* x[i,j] for i in 1:n, j in i+1:n) <= z ))
    println("==> ajout de la contrainte SP1")

    # On relance la recherche de solution de probleme maitre
    optimize!(m)
  
    z_etoile = objective_value(m)
    #print("==> nouvelle z* =",z_etoile,"<=======")
    x_etoile=value.(x)
    y_etoile=value.(y)
    # On met à jour les objectifs des sous-problèmes
    @objective(m1, Max, sum((l[i,j]+sigma1[i,j]*(lh[i]+lh[j])) * x_etoile[i,j] for i in 1:n, j in i+1:n))
    #nb+=1
    #println("========>nous sommes rentrés dans la boucle:",nb,"<=========")
    # On relance les optimisations
    optimize!(m1)
    z1=objective_value(m1) # Nouvelles valeurs des valeurs optimales des sous-problèmes
    vsigma1=value.(sigma1)

end
y_etoile=value.(y)
println("la valeur est = ",z_etoile)