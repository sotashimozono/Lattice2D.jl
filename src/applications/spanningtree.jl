"""
    Spanning Trees Module

This module implements uniform random spanning tree generation using Wilson's algorithm.
A spanning tree is a subgraph that connects all vertices without cycles.
"""

"""
    SpanningTree

Represents a spanning tree on a lattice.

# Fields
- `lattice::Lattice`: The underlying lattice
- `tree_edges::Set{Tuple{Int,Int}}`: Edges in the spanning tree
- `parent::Vector{Int}`: Parent of each vertex in the tree (-1 for root)
- `root::Int`: Root vertex of the tree
"""
struct SpanningTree{T,B,I}
    lattice::Lattice{T,Float64,B,I}
    tree_edges::Set{Tuple{Int,Int}}
    parent::Vector{Int}
    root::Int
end

export SpanningTree

"""
    wilson_algorithm(lattice::Lattice, root::Int=1; rng=Random.default_rng())

Generate a uniform random spanning tree using Wilson's algorithm.

Wilson's algorithm performs loop-erased random walks from each vertex
to build a uniform random spanning tree.

# Arguments
- `lattice::Lattice`: The lattice
- `root::Int`: Root vertex (default: 1)
- `rng`: Random number generator (optional)

# Returns
- `SpanningTree`: The generated spanning tree

"""
function wilson_algorithm(lattice::Lattice{T,F,B,I}, root::Int=1; 
                         rng=Random.default_rng()) where {T,F,B,I}
    @assert 1 <= root <= lattice.N "Root must be within lattice bounds"
    
    N = lattice.N
    in_tree = falses(N)
    parent = fill(-1, N)
    tree_edges = Set{Tuple{Int,Int}}()
    
    # Start with root
    in_tree[root] = true
    parent[root] = -1
    
    # Add vertices one by one
    for start in 1:N
        if in_tree[start]
            continue
        end
        
        # Perform loop-erased random walk from start
        current = start
        path = Dict{Int,Int}()  # maps vertex to next vertex in path
        
        while !in_tree[current]
            neighbors = lattice.nearest_neighbors[current]
            next = rand(rng, neighbors)
            path[current] = next
            current = next
        end
        
        # Add the loop-erased path to tree
        current = start
        while !in_tree[current]
            next = path[current]
            in_tree[current] = true
            parent[current] = next
            
            # Add edge (ensuring consistent ordering)
            edge = current < next ? (current, next) : (next, current)
            push!(tree_edges, edge)
            
            current = next
        end
    end
    
    return SpanningTree{T,B,I}(lattice, tree_edges, parent, root)
end
export wilson_algorithm

"""
    is_connected(tree::SpanningTree)

Check if the spanning tree connects all vertices.

# Returns
- `Bool`: true if tree connects all vertices

"""
function is_connected(tree::SpanningTree)
    N = tree.lattice.N
    visited = falses(N)
    
    # BFS from root
    queue = [tree.root]
    visited[tree.root] = true
    
    while !isempty(queue)
        v = popfirst!(queue)
        
        # Find children
        for i in 1:N
            if tree.parent[i] == v && !visited[i]
                visited[i] = true
                push!(queue, i)
            end
        end
    end
    
    return all(visited)
end
export is_connected

"""
    tree_depth(tree::SpanningTree, vertex::Int)

Calculate the depth (distance from root) of a vertex.

# Arguments
- `tree::SpanningTree`: The spanning tree
- `vertex::Int`: Vertex index

# Returns
- `Int`: Depth of vertex

"""
function tree_depth(tree::SpanningTree, vertex::Int)
    @assert 1 <= vertex <= tree.lattice.N "Vertex must be within lattice bounds"
    
    depth = 0
    current = vertex
    
    while tree.parent[current] != -1
        current = tree.parent[current]
        depth += 1
    end
    
    return depth
end
export tree_depth

"""
    tree_height(tree::SpanningTree)

Calculate the maximum depth (height) of the tree.

# Returns
- `Int`: Height of the tree

"""
function tree_height(tree::SpanningTree)
    return maximum(tree_depth(tree, v) for v in 1:tree.lattice.N)
end
export tree_height

"""
    path_to_root(tree::SpanningTree, vertex::Int)

Find the path from a vertex to the root.

# Arguments
- `tree::SpanningTree`: The spanning tree
- `vertex::Int`: Starting vertex

# Returns
- `Vector{Int}`: Path from vertex to root

"""
function path_to_root(tree::SpanningTree, vertex::Int)
    @assert 1 <= vertex <= tree.lattice.N "Vertex must be within lattice bounds"
    
    path = [vertex]
    current = vertex
    
    while tree.parent[current] != -1
        current = tree.parent[current]
        push!(path, current)
    end
    
    return path
end
export path_to_root

"""
    tree_distance(tree::SpanningTree, u::Int, v::Int)

Calculate the distance between two vertices in the tree.

# Arguments
- `tree::SpanningTree`: The spanning tree
- `u::Int`: First vertex
- `v::Int`: Second vertex

# Returns
- `Int`: Tree distance (number of edges)

"""
function tree_distance(tree::SpanningTree, u::Int, v::Int)
    path_u = path_to_root(tree, u)
    path_v = path_to_root(tree, v)
    
    # Find lowest common ancestor
    set_u = Set(path_u)
    lca = -1
    for vertex in path_v
        if vertex in set_u
            lca = vertex
            break
        end
    end
    
    # Distance = dist(u, lca) + dist(v, lca)
    dist_u = findfirst(x -> x == lca, path_u) - 1
    dist_v = findfirst(x -> x == lca, path_v) - 1
    
    return dist_u + dist_v
end
export tree_distance

"""
    subtree_size(tree::SpanningTree, vertex::Int)

Calculate the size of the subtree rooted at a vertex.

# Arguments
- `tree::SpanningTree`: The spanning tree
- `vertex::Int`: Root of subtree

# Returns
- `Int`: Number of vertices in subtree

"""
function subtree_size(tree::SpanningTree, vertex::Int)
    @assert 1 <= vertex <= tree.lattice.N "Vertex must be within lattice bounds"
    
    # Count descendants
    count = 1  # count self
    
    for i in 1:tree.lattice.N
        if tree.parent[i] == vertex
            count += subtree_size(tree, i)
        end
    end
    
    return count
end
export subtree_size

"""
    leaf_vertices(tree::SpanningTree)

Find all leaf vertices (vertices with no children).

# Returns
- `Vector{Int}`: List of leaf vertex indices

"""
function leaf_vertices(tree::SpanningTree)
    N = tree.lattice.N
    has_child = falses(N)
    
    for i in 1:N
        p = tree.parent[i]
        if p != -1
            has_child[p] = true
        end
    end
    
    leaves = Int[]
    for i in 1:N
        if !has_child[i]
            push!(leaves, i)
        end
    end
    
    return leaves
end
export leaf_vertices

"""
    average_path_length(tree::SpanningTree; n_samples::Int=100, rng=Random.default_rng())

Estimate the average path length between random vertex pairs.

# Arguments
- `tree::SpanningTree`: The spanning tree
- `n_samples::Int`: Number of random pairs to sample (default: 100)
- `rng`: Random number generator (optional)

# Returns
- `Float64`: Average path length

"""
function average_path_length(tree::SpanningTree; n_samples::Int=100, rng=Random.default_rng())
    N = tree.lattice.N
    total_dist = 0
    
    for _ in 1:n_samples
        u = rand(rng, 1:N)
        v = rand(rng, 1:N)
        total_dist += tree_distance(tree, u, v)
    end
    
    return total_dist / n_samples
end
export average_path_length

"""
    degree_sequence(tree::SpanningTree)

Calculate the degree sequence of vertices in the tree.

# Returns
- `Vector{Int}`: Degree of each vertex

"""
function degree_sequence(tree::SpanningTree)
    N = tree.lattice.N
    degrees = zeros(Int, N)
    
    for edge in tree.tree_edges
        u, v = edge
        degrees[u] += 1
        degrees[v] += 1
    end
    
    return degrees
end
export degree_sequence

"""
    _prune_leaves!(tree::SpanningTree)

Internal function: Remove all leaf vertices from the tree (modifies tree in place).

# Returns
- `Int`: Number of leaves removed

**Note**: This is an internal function and may change in future versions.
"""
function _prune_leaves!(tree::SpanningTree)
    leaves = leaf_vertices(tree)
    
    for leaf in leaves
        if leaf == tree.root
            continue  # Don't remove root
        end
        
        # Remove edge to parent
        p = tree.parent[leaf]
        edge = leaf < p ? (leaf, p) : (p, leaf)
        delete!(tree.tree_edges, edge)
        tree.parent[leaf] = -1
    end
    
    return length(leaves)
end

"""
    generate_maze(lattice::Lattice; rng=Random.default_rng())

Generate a maze using a spanning tree. 
The maze is a spanning tree where passages correspond to tree edges.

# Arguments
- `lattice::Lattice`: The lattice (should have OBC for a proper maze)
- `rng`: Random number generator (optional)

# Returns
- `SpanningTree`: The maze as a spanning tree

"""
function generate_maze(lattice::Lattice; rng=Random.default_rng())
    # Choose corner as root for better maze structure
    root = 1
    return wilson_algorithm(lattice, root, rng=rng)
end
export generate_maze

"""
    critical_path(tree::SpanningTree)

Find the longest path in the tree (diameter).

# Returns
- `Tuple{Vector{Int}, Int}`: (path, length)

"""
function critical_path(tree::SpanningTree)
    N = tree.lattice.N
    
    # First BFS from root to find farthest vertex
    depths = [tree_depth(tree, v) for v in 1:N]
    v1 = argmax(depths)
    
    # Reconstruct tree with v1 as root
    # This is a simplified approach - just find farthest from v1
    max_dist = 0
    v2 = v1
    
    for v in 1:N
        dist = tree_distance(tree, v1, v)
        if dist > max_dist
            max_dist = dist
            v2 = v
        end
    end
    
    # Build path from v1 to v2
    path_v1 = path_to_root(tree, v1)
    path_v2 = path_to_root(tree, v2)
    
    # Find LCA and construct path
    set_v1 = Set(path_v1)
    lca_idx = findfirst(v -> v in set_v1, path_v2)
    
    if lca_idx === nothing
        # This should never happen in a connected tree
        @warn "Could not find LCA in critical_path - tree may be disconnected"
        return Int[], 0
    end
    
    lca = path_v2[lca_idx]
    
    # Path from v1 to lca
    lca_idx_v1 = findfirst(v -> v == lca, path_v1)
    path1 = reverse(path_v1[1:lca_idx_v1])
    
    # Path from lca to v2
    path2 = path_v2[1:lca_idx]
    
    # Combine (excluding duplicate lca)
    full_path = vcat(path1[1:end-1], path2)
    
    return full_path, length(full_path) - 1
end
export critical_path
