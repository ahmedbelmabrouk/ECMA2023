using SimpleWeightedGraphs
using IterTools
using Munkres
using LinearAlgebra
using SparseArrays
using Graphs
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
function matching_cost(indexes,graph)
    return sum(graph[j[1],j[2]] for j in indexes)
end
function transform(list_indexes)
    return[(i,list_indexes[i]) for i in 1:length(list_indexes)]
end
function create_multigraph(tsp,mst,matching_indexes,odd_degree)
    multigraph=Multigraph(SimpleWeightedGraphs.nv(tsp))
    for edge in collect(SimpleWeightedGraphs.edges(mst))
        Multigraph.add_edge!(multigraph,edge)
    end
    for pair in matching_indexes
        Multigraph.add_edge!(multigraph,pair[1],pair[2])
    end
    return multigraph
end
