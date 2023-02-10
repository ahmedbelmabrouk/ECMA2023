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
l=distance(n,coordinates)
donnee=n,L,B,K,W_v,w_v,W,coordinates,l
function U1(donnee,x)
    n,L,B,K,W_v,w_v,W,coordinates,l=donnee
    m1 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_MIPDISPLAY" =>0))
    @variable(m1, sigma1[i in 1:n,j in 1:n] >=0.01)
    @constraint(m1, sum(sigma1[i,j] for i in 1:n,j in i+1:n) <= L)
    @constraint(m1 , [i in 1:n, j in 1:n] , sigma1[i,j] <= 3)
    
    
    
    @objective(m1,Max, sum((l[i,j]+sigma1[i,j]*(lh[i]+lh[j]))*x[i,j] for i=1:n, j=i+1:n))
    optimize!(m1)
    l1=zeros(n,n)
    for i in 1:n
        for j=1:n
            l1[i,j]=l[i,j]+value.(sigma1)[i,j]*(lh[i]+lh[j])
        end
    end
    return l1
end
function U2(donnee)
    n,L,B,K,W_v,w_v,W,coordinates,l=donnee
    m2 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_MIPDISPLAY" =>0))
    @variable(m2, sigma2[i in 1:n] >=0.05)
    @constraint(m2, (sum(sigma2[i] for i in 1:n) <= W))
    @constraint(m2, [i in 1:n] , sigma2[i] <= W_v[i] )
    
    optimize!(m2)
    w_v2=[w_v[i]*(1+value.(sigma2)[i]) for i in 1:n]
    return w_v2
end
function solution_triviale(donnee)

    n,L,B,K,W_v,w_v,W,coordinates,l=donnee
    
    w_v2=U2(donnee)
    m=Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_MIPDISPLAY" =>0))
    @variable(m,y[i=1:n,k=1:K],Bin)

    
    #Définition des contraintes

    @constraint(m,[k in 1:K], sum(w_v2[i] * y[i,k] for i in 1:n ) <= B)
    @constraint(m,[i in 1:n],sum(y[i,k] for k in 1:K ) == 1)
    optimize!(m)

    return value.(y)
end
function coexistance(solution,donnee)
    n=donnee[1]
    x=Matrix{Int64}(zeros(n,n))
    for i in 1:n
        for j in 1:n 
            coex=[solution[i,k]==solution[j,k] for k in 1:donnee[4]]
            if coex==[true for k in 1:donnee[4]]
                x[i,j]=1
            else 
                x[i,j]=0
            end
        end
    end
    return x
end
function objective(solution,donnee)
    x=coexistance(solution,donnee)
    l1=U1(donnee,x)
    n=donnee[1]
    v=0
    for i in 1:n
        for j in 1:n
            v+=l1[i,j]*x[i,j]
        end
    end
    return v
end

function cluster(solution, i ,donnee)
    K=donnee[4]
    for j=1:K
        if solution[i,j]==1
            return j 
        end
    end
end
function valid(solution, donnee)
    n=donnee[1]
    B=donnee[3]
    w_v=donnee[6]
    K=donnee[4]
    verif1=[sum(w_v[i]*solution[i,k] for i in 1:n)<=B for k=1:K ]
    verif2=[sum(solution[i,k] for k in 1:K)==1 for i in 1:n]
    verif=[verif1;verif2]
    for i in verif 
        if i==false
            return false
        end
    end
    return true
end


function permut(solution,i1,i2,donnee)
    k1=cluster(solution,i1,donnee)
    k2=cluster(solution,i2,donnee)
    result=copy(solution)
    result[i1,k1]=0
    result[i1,k2]=1
    result[i2,k2]=0
    result[i2,k1]=1
    return result
end
function ls_permut(donnee,nb_duration_max=1800,nb_iter_max=10000)
    iter=0
    nb_cons_reject=0
    nb_move=0
    start=time()
    duration=time()-start
    finished=(iter>nb_iter_max)||(nb_cons_reject>10000)||(duration>nb_duration_max)
    cursol=solution_triviale(donnee)
    
    bestsol=cursol
    n=donnee[1]


    while !finished 
        iter+=1
        i1=rand(1:n)
        i2=rand([i for i in 1:n if i!=i1])
        testsol=permut(cursol,i1,i2,donnee)
        while !valid(testsol,donnee)
            i1=rand(1:n)
            i2=rand([i for i in 1:n if i!=i1])
            testsol=permut(cursol,i1,i2,donnee)
        end
        println(objective(testsol,donnee))
        if objective(testsol,donnee)<objective(cursol,donnee)
            nb_cons_reject=0
            nb_move+=1
            copy!(cursol,testsol)
            if objective(cursol,donnee)<objective(bestsol,donnee)
                copy!(bestsol,cursol)
                println("bestsol====$bestsol")
                println("cost=====$(objective(bestsol,donnee))")

            end
        else
            nb_cons_reject+=1
        end
        duration=time()-start
        finished=(iter>nb_iter_max)||(duration>nb_duration_max)
        #finished=(duration>nb_duration_max)
    end
    return [objective(bestsol,donnee),iter,nb_cons_reject,duration,nb_move]
end
data=["data/10_ulysses_6.tsp","data/10_ulysses_9.tsp","data/14_burma_3.tsp","data/14_burma_6.tsp","data/14_burma_9.tsp","data/22_ulysses_3.tsp","data/22_ulysses_6.tsp","data/22_ulysses_9.tsp","data/26_eil_3.tsp"]
for path in data 
    open("results_descent.txt","a") do file 
        include(path)
        l=distance(n,coordinates)
        donnee=n,L,B,K,W_v,w_v,W,coordinates,l
        ms=ls_permut(donnee,300)
        
        ch=append!([string(n),string(K)],[string(i) for i in ms])
        println(file,ch)
    end
end
include(data[1])
        l=distance(n,coordinates)
        donnee=n,L,B,K,W_v,w_v,W,coordinates,l
        ms=ls_permut(donnee,900,1000)
        
        ch=[string(i) for i in ms]
        println(ch)
