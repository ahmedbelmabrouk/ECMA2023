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


#***************Plan_Coupants****************************
function PlansCP(n,L,B,K,W_v,w_v,W,coordinates,l)
#**********Probleme maitre************
m = Model(CPLEX.Optimizer)

@variable(m, x[i in 1:n, j in 1:n] , Bin)
@variable(m,y[i=1:n,k=1:K],Bin)
@variable(m, z >=0)

#Définition des contraintes
@constraint(m, (sum(l[i,j] * x[i,j] for i in 1:n , j in 1:n )) <= z )
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K ],x[i,j] + y[i,k] - y[j,k] <= 1) #lien entre x et y 1/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K ],x[i,j] - y[i,k] + y[j,k] <= 1) #lien entre x et y 2/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K ],-x[i,j] + y[i,k] + y[j,k] <= 1) #lien entre x et y 3/3
@constraint(m,[k in 1:K], sum(w_v[i] * y[i,k] for i in 1:n ) <= B)
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
@constraint(m1, sum(sigma1[i,j] for i in 1:n,j in 1:n) <= L)
@constraint(m1 , [i in 1:n, j in 1:n] , sigma1[i,j] <= 3)
@objective(m1, Max, sum((l[i,j]+sigma1[i,j]*(lh[i]+lh[j])) * x_etoile[i,j] for i in 1:n, j in 1:n))
optimize!(m1)
z1=objective_value(m1)
vsigma1=value.(sigma1)


#***********k sous problme 2************
function SP2(k,y_etoile)
    m2k = Model(CPLEX.Optimizer)
    @variable(m2k, sigma2k[i in 1:n] >=0)
    @constraint(m2k, (sum(sigma2k[i] for i in 1:n) <= W))
    @constraint(m2k, [i in 1:n] , sigma2k[i] <= W_v[i] )
    @objective(m2k, Max, sum((w_v[i] *(1+ sigma2k[i])* y_etoile[i,k]) for i in 1:n))
    optimize!(m2k)
    z2k=objective_value(m2k)
    vsigma2k=value.(sigma2k)
    return(z2k,vsigma2k)
end
z2=[]
vsigma2=[]
for k in 1:K
    push!(z2,SP2(k,y_etoile)[1])
    push!(vsigma2,SP2(k,y_etoile)[2])
end
"""
m2=[ Model(CPLEX.Optimizer) for k in 1:K]
z2=[0.1 for k in 1:K]
vsigma2=[[] for k in 1:K]
for k in 1:K

m2[k] = Model(CPLEX.Optimizer)

@variable(m2[k], sigma2[i in 1:n,k in 1:K] >=0)
@constraint(m2[k], (sum(sigma2[i ,k] for i in 1:n) <= W))
@constraint(m2[k], [i in 1:n] , sigma2[i,k] <= W_v[i] )
@objective(m2[k], Max, sum((w_v[i] *(1+ sigma2[i,k])*y_etoile[i,k]) for i in 1:n))
optimize!(m2[k])
z2[k]=objective_value(m2[k])
vsigma2[k]=value.(sigma2[:,k])
end
"""
nb=0
while z1 > z_etoile || maximum(z2) > B 
    if z1 > z_etoile > 1e-4
        @constraint(m, (sum((l[i,j]+ vsigma1[i,j]*(lh[i]+lh[j]))* x[i,j] for i in 1:n, j in 1:n) <= z ))
        println("==> ajout de la contrainte SP1")
    end
    for k in 1:K
        if z2[k] > B
            @constraint(m,sum(w_v[i]* (1+ vsigma2[k][i])* y[i,k] for i in 1:n, j in 1:n)<=B)
            println("==> ajout de keme contraintes SP2")
        end
    end
    
    
    # On relance la recherche de solution de probleme maitre
    optimize!(m)
  
    z_etoile = objective_value(m)
    print("==> nouvelle z* =",z_etoile,"<=======")
    x_etoile=value.(x)
    y_etoile=value.(y)
    # On met à jour les objectifs des sous-problèmes
    @objective(m1, Max, sum((l[i,j]+sigma1[i,j]*(lh[i]+lh[j])) * x_etoile[i,j] for i in 1:n, j in 1:n))
    nb+=1
    println("========>nous sommes rentrés dans la boucle:",nb,"<=========")
    # On relance les optimisations
    optimize!(m1)
    z1=objective_value(m1) # Nouvelles valeurs des valeurs optimales des sous-problèmes
    vsigma1=value.(sigma1)
    z2=[]
    vsigma2=[]
    for k in 1:K
        push!(z2,SP2(k,y_etoile)[1])
        push!(vsigma2,SP2(k,y_etoile)[2])
    end
    """
    for k in 1:K
        @objective(m2[k], Max, sum((w_v[i] *(1+ sigma2[i,k])* y_etoile[i,k]) for i in 1:n))
        optimize!(m2[k])
        vsigma2[k]=value.(sigma2[:,k])
        z2[k]=objective_value(m2[k]) # Nouvelles valeurs des valeurs optimales des sous-problèmes
    end
    """
end

return (z_etoile)
end
val=PlansCP(n,L,B,K,W_v,w_v,W,coordinates,l)
println("la valeur est = ",val)
"""
#******************MAIN****************
#liste des noms des instances à tester
L=[]
for i in 20:20:200
    for v in ["BAY","COL","NY"]
        ch=string(i) * "_USA-road-d." * v * ".gr"
        append!(L,[ch])
    end
end
for i in 200:50:1000
    for v in ["NY","BAY","COL"]
        ch=string(i) * "_USA-road-d." * v * ".gr"
        append!(L,[ch])
    end
end
for i in 1000:100:1200
    for v in ["NY","BAY","COL"]
        ch=string(i) * "_USA-road-d." * v * ".gr"
        append!(L,[ch])
    end
end
#***************************
tab_res3=[]
for nom in L[1:30]

    (n,s,t,S,d1,d2,p,ph,Mat)=LireInstance(nom)
    
    res= PlansCoupants((n,s,t,S,d1,d2,p,ph,Mat))
    t= @elapsed PlansCoupants((n,s,t,S,d1,d2,p,ph,Mat))
    append!(tab_res3,[[res,t]])
end
    
print(tab_res3)

""""