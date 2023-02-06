using SimpleWeightedGraphs
using IterTools
using Munkres
using LinearAlgebra
using SparseArrays
using Graphs
using Hungarian
using Multigraphs
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


function creating_tsp(donnee)
    n=donnee[1]
    l=donnee[9]
    TSP=SimpleWeightedGraph(n+1)
    for i=1:n
        for j=1:n
            SimpleWeightedGraphs.add_edge!(TSP,i,j,l[i,j])
        end
        SimpleWeightedGraphs.add_edge!(TSP,i,n+1,10000)
    end
    return TSP
end

function minimum_spanning_tree(tsp,donnee) # calcule l'arbre couvrant de poids minimum. 
	tree=SimpleWeightedGraphs.kruskal_mst(tsp)
    mst=SimpleWeightedGraph(length(tree)+1)
    for i=1:length(tree)
        SimpleWeightedGraphs.add_edge!(mst,tree[i])  
    end 
    return mst
end
function odd_degree(mst,n)
    return [i for i in SimpleWeightedGraphs.vertices(mst) if SimpleWeightedGraphs.degree(mst,i)%2==1]
end
function bipartite_set(odd_degree_vertices)
    tt=round(Int64,length(odd_degree_vertices)/2)
    h=collect(subsets(odd_degree_vertices,tt))
    return collect([i for i in collect(subsets(odd_degree_vertices,tt)) ])
end
function graph_biparti(tsp,bipartite_sett,odd_degree_vertices)
    bipartite_graphs=[]
    vertices_sets=[]
    for vertex_set1 in bipartite_sett
        vertex_set1=sort!(vertex_set1)
        vertex_set2=[]
        for vertex in odd_degree_vertices
            if !(vertex in vertex_set1)
                push!(vertex_set2,vertex)
            end
        end
        vertex_set2=Vector{Int64}(vertex_set2)
        matrix=[10000.00 for _ in 1:length(vertex_set2), _ in 1:length(vertex_set1)]
        for i=1:length(vertex_set1)
            for j=1:length(vertex_set2)
                if vertex_set1[i] < vertex_set2[j]
                    matrix[i,j]=SimpleWeightedGraphs.weights(tsp)[vertex_set1[i],vertex_set2[j]]
                else
                    matrix[i,j]=SimpleWeightedGraphs.weights(tsp)[vertex_set2[j],vertex_set1[i]]
                end
            end
        end
        push!(bipartite_graphs,matrix)
        push!(vertices_sets,[vertex_set1,vertex_set2])
    end
    return (bipartite_graphs,vertices_sets)
end
function transform(list_indexes)
    return[(i,list_indexes[i]) for i in 1:length(list_indexes)]
end
function matching_cost(indexes,graph)
    return sum(graph[j[1],j[2]] for j in indexes)
end
function min_matching(tsp,graphs)
    n=SimpleWeightedGraphs.size(tsp)[1]
    minii=10^12
    min_index=1
    min_matching_indexes=[]
    for (index,graph) in  collect(enumerate(graphs[1]))
        
        matching_indexes=transform(hungarian(graph)[1])
        cost=matching_cost(matching_indexes,graph)
        if cost<minii
            minii=cost
            min_index=index
            min_matching_indexes=matching_indexes
        end
    end
    println(min_index,min_matching_indexes)
    matching_index=[[] for _ in 1:length(min_matching_indexes)]
    for (index,vertex_set) in collect(enumerate(min_matching_indexes))
        append!(matching_index[index],graphs[2][min_index][1][vertex_set[1]])
        append!(matching_index[index],graphs[2][min_index][2][vertex_set[2]])
    end
    return matching_index

end


function create_multigraph(tsp,mst,matching_indexes,odd_degree,donnee)
    nn=SimpleWeightedGraphs.nv(tsp)
    multigraph=SimpleGraph(nn)
    for edge in collect(SimpleWeightedGraphs.edges(mst))
        add_edge!(multigraph,src(edge),dst(edge))
        println(edge)
    end
    for pair in matching_indexes
        add_edge!(multigraph,pair[1],pair[2])
        println(pair)
    end
    return multigraph
end
function weight(tour,n,w_v)
    w=0
    for i in tour
        if i!=n+1
            w+=w_v[i]
        end
    end
    return w 
end

include("data/10_ulysses_9.tsp")
l=distance(n,coordinates)
donnee=n,L,B,K,W_v,w_v,W,coordinates,l

tsp=creating_tsp(donnee)
mst=minimum_spanning_tree(tsp,donnee)
oo=odd_degree(mst,n)
ss=bipartite_set(oo)
graphs=graph_biparti(tsp,ss,oo)
matching_indexes=min_matching(tsp,graphs)
mg=create_multigraph(tsp,mst,matching_indexes,oo,donnee)

data=["data/10_ulysses_3.tsp","data/10_ulysses_6.tsp","data/10_ulysses_9.tsp","data/14_burma_3.tsp","data/14_burma_6.tsp","data/14_burma_9.tsp","data/22_ulysses_3.tsp","data/22_ulysses_6.tsp","data/22_ulysses_9.tsp","data/26_eil_3.tsp"]
for path in data 
    include(path)
    l=distance(n,coordinates)
    donnee=n,L,B,K,W_v,w_v,W,coordinates,l
    tsp=creating_tsp(donnee)
    mst=minimum_spanning_tree(tsp,donnee)
    oo=odd_degree(mst,n)
    ss=bipartite_set(oo)
    graphs=graph_biparti(tsp,ss,oo)
    matching_indexes=min_matching(tsp,graphs)
    mg=create_multigraph(tsp,mst,matching_indexes,oo,donnee)
    
    
end
include("data/10_ulysses_3.tsp")


    B=donnee[3]
    w_v=donnee[6]
    append!(w_v,0)
    n=donnee[1]
    weight_tour=0
    
    tour=[]
    multigraph=copy(mg)
    temp_graph=Graph(0)
    graph_vertices=vertices(mg)
    current_vertex=graph_vertices[1]
    append!(tour,current_vertex)
    weight_tour=weight(tour,n,w_v)
    while ne(multigraph)!=0
        
        for edge in [edge for edge in collect(edges(multigraph)) if src(edge)==current_vertex ||dst(edge)==current_vertex]
            temp_graph=deepcopy(multigraph)
            println("edge teste=$edge")
            rem_edge!(temp_graph,src(edge),dst(edge))
            if weight_tour+w_v[dst(edge)]>B
                println(tour)
                break
            end
            println(collect(edges(temp_graph)))
            if is_connected(temp_graph)
                append!(tour,dst(edge))
                weight_tour=weight(tour,n,w_v)
                current_vertex=dst(edge)
                rem_edge!(multigraph,src(edge),dst(edge))
                println("tour=====$tour")
                println("==="^70)
                break
            else
                append!(tour,dst(edge))
                weight_tour=weight(tour,n,w_v)
                current_vertex=dst(edge)
                rem_edge!(multigraph,src(edge),dst(edge))
                println("tour=====$tour")
                println("==="^70)
                isolates=[v for v in vertices(multigraph) if length(all_neighbors(multigraph,v))==0]
                for v in isolates 
                    rem_vertex!(multigraph,v)
                end
                
            end
            
        end
    end


    for u in vertices(g)
        label[u] != 0 && continue
        label[u] = u
        Q = Vector{Int64}()
        push!(Q, u)
        while !isempty(Q)
            src = popfirst!(Q)
            for vertex in all_neighbors(g, src)
                if label[vertex] == 0
                    push!(Q, vertex)
                    label[vertex] = u
                end
            end
        end
    end

