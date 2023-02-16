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
MOI.set(m, MOI.NumberOfThreads(), 1)
@variable(m, x[i in 1:n, j in 1:n] , Bin)
@variable(m,y[i=1:n,k=1:K],Bin)
@variable(m, z >=0)

#Définition des contraintes
@constraint(m, (sum(l[i,j] * x[i,j] for i in 1:n , j in i+1:n )) <= z )
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j ],x[i,j] + y[i,k] - y[j,k] <= 1) #lien entre x et y 1/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j ],x[i,j] - y[i,k] + y[j,k] <= 1) #lien entre x et y 2/3
@constraint(m,[i in 1:n ,j in 1:n , k in 1:K;i!=j ],-x[i,j] + y[i,k] + y[j,k] <= 1) #lien entre x et y 3/3
@constraint(m,[k in 1:K], sum(w_v[i] * y[i,k] for i in 1:n ) <= B)
@constraint(m,[i in 1:n],sum(y[i,k] for k in 1:K ) == 1)

#Définition de l'objectif
@objective(m, Min, z)

#**********sous probleme 1************
m1 = Model(CPLEX.Optimizer)

@variable(m1, sigma1[i in 1:n,j in 1:n] >=0)
@constraint(m1, sum(sigma1[i,j] for i in 1:n,j in i+1:n) <= L)
@constraint(m1 , [i in 1:n, j in 1:n] , sigma1[i,j] <= 3)

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
function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end

function my_cb_function(cb_data::CPLEX.CallbackContext, context_id::Clong)
    if isIntegerPoint(cb_data, context_id)
        CPLEX.load_callback_variable_primal(cb_data, context_id) #Cette ligne doit être appelée avant de pouvoir récupérer la solution entière ayant entraîné l’appel du callback
        # On récupère la valeur de z
        
        x_etoile = callback_value.(cb_data, x)
        y_etoile = callback_value.(cb_data, y)
        z_val = callback_value(cb_data, z)
        #println("===========> z*=",z_val)
        
        @objective(m1, Max, sum((l[i,j]+sigma1[i,j]*(lh[i]+lh[j])) * x_etoile[i,j] for i in 1:n, j in i+1:n))
        optimize!(m1)
        z1=objective_value(m1)
        vsigma1=value.(sigma1)
        
        z2=[]
        vsigma2=[]
        for k in 1:K
            push!(z2,SP2(k,y_etoile)[1])
            #println("!<!<!<!<!<",z2[k])
            push!(vsigma2,SP2(k,y_etoile)[2])
            #println("!<!<!<!<!<",vsigma2[k])
        end
        #println(y_etoile)
        #println(vsigma2)

        if (z1 - z_val) > 1e-4
            cstr1 = @build_constraint (sum((l[i,j]+ vsigma1[i,j]*(lh[i]+lh[j]))* x[i,j] for i in 1:n, j in i+1:n) <= z )
            MOI.submit(m, MOI.LazyConstraint(cb_data), cstr1)
            #println("==> ajout de la contrainte SP1")
        end
        if maximum(z2)-B > 1e-4
            for k in 1:K
                #println(">>>>>>>>>>>VAL VIOL:",z2[k] - B)
                if z2[k] - B > 1e-4
                    cstr2 = @build_constraint (sum(w_v[i]* (1+ vsigma2[k][i])* y[i,k] for i in 1:n)<=B)
                    MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2)
                    #println("==> ajout de keme contraintes SP2",cstr2)
                    
                end
            end
        end
    end
end

# On précise que le modèle doit utiliser notre fonction de callback
MOI.set(m, CPLEX.CallbackFunction(), my_cb_function)
optimize!(m)
z_etoile = objective_value(m)
println("Objective value: ", z_etoile)
#return (z_etoile)
