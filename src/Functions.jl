using StatsBase  # for countmap or not
using Statistics
using Base.Threads


# using Threads
## Input: input dataset in hMetis format
# Outout: hypergraph array
function ReadInp(input)

	io = open(input, "r")

	ar = Any[]
	while !eof(io)
		rr = zeros(Int, 0)
		ln = readline(io)
		sp = split(ln)

		for kk ∈ 1:length(sp)
			r = parse(Int, sp[kk])
			append!(rr, r)
		end #kk
		push!(ar, rr)

	end

	ar = deleteat!(ar, 1)

	return ar

end #end of function

## Input: hypergraph array
# Outout: sparse incidence matrix
function INC(ar)

	col = zeros(Int, 0)
	row = zeros(Int, 0)



	for iter ∈ 1:length(ar)
		cc = (iter) * ones(Int, length(ar[iter]))
		rr = ar[iter]

		append!(col, cc)
		append!(row, rr)
	end

	row = row

	val = ones(Float64, length(row))

	mat = sparse(col, row, val)

	return mat
end

## Input: hypergraph array
# Output: number of nodes in hypergraph
function mxF(ar)

	mx2 = Int(0)
	aa = Int(0)

	for i ∈ 1:length(ar)

		mx2 = max(aa, maximum(ar[i]))
		aa = mx2

	end
	return mx2

end

## Input: hypergraph array, hyperedge weights
# Output: sparse simple graph
function StarW(ar, W)

	mx = mxF(ar)

	sz = length(ar)
	col = zeros(Int32, 0)
	val = zeros(Float32, 0)
	row = zeros(Int32, 0)

	for iter ∈ 1:length(ar)
		LN = length(ar[iter])
		cc = (iter + mx) * ones(Int, LN)
		vv = (W[iter] / LN) * ones(Int, LN)

		rr = ar[iter]
		append!(col, cc)

		append!(row, rr)

		append!(val, vv)
	end

	mat = sparse(row, col, val, mx + sz, mx + sz)

	A = mat + mat'

	return A

end

function CliqueW(ar, W)
	# 1) Determine maximum node index
	mx = mxF(ar)

	# 2) We'll build a node×node adjacency matrix in dimension (mx × mx)
	#    row, col, val for the final sparse matrix
	row = Int[]
	col = Int[]
	val = Float64[]

	# 3) For each hyperedge e, connect all pairs (u, v) of nodes in that hyperedge
	for (e_idx, edge_nodes) in enumerate(ar)
		LN = length(edge_nodes)
		if LN < 2
			continue
		end
		# Distribute the hyperedge weight among all pairs
		pair_weight = W[e_idx] / (LN * (LN - 1) / 2)

		# For each unordered pair (i, j)
		for i in 1:(LN-1)
			for j in (i+1):LN
				n1 = edge_nodes[i]
				n2 = edge_nodes[j]
				push!(row, n1)
				push!(col, n2)
				push!(val, pair_weight)
				# Symmetrize
				push!(row, n2)
				push!(col, n1)
				push!(val, pair_weight)
			end
		end
	end

	# 4) Build the sparse adjacency matrix
	A_clique = sparse(row, col, val, mx, mx)

	return A_clique
end

function CliqueSampler(hyperedge::Vector{Int}, edge_weights::Float64; sample_num::Int = 1)

	len = length(hyperedge)
	I = Int[]
	J = Int[]
	V = Float64[]

	if len < 2
		# raise an error
		error("The number of node in an hyperedge must be at least 2.")
	else
		# Sort the hyperedge node by node weights
		equivalent_node_weights = 2 * edge_weights / (len - 1)
		node_weights = fill(equivalent_node_weights, len)
		cumulative_sum = cumsum(node_weights)
		total_weights = cumulative_sum[end]

		if len == 2
			push!(I, hyperedge[1])
			push!(J, hyperedge[2])
			push!(V, node_weights[1] * node_weights[2] / total_weights)
		else
			o = Base.Order.ord(isless, identity, false, Base.Order.Forward) # for ascending order sorting

			for _ in 1:sample_num

				for i in 1:(len-1)
					new_edeg_weights = node_weights[i] * (total_weights - cumulative_sum[i]) / total_weights / sample_num

					# sample the target node
					r = rand() * (total_weights - cumulative_sum[i]) + cumulative_sum[i]
					j = searchsortedfirst(cumulative_sum, r, 1, len, o)

					i = hyperedge[i]
					j = hyperedge[j]

					if i > j
						temp = j
						j = i
						i = temp
					end


					push!(I, i)
					push!(J, j)
					push!(V, new_edeg_weights)

				end # for

			end # for


		end # if
	end # if

	return I, J, V

end
function CliqueExpansion(hyperedge, edge_weights)

	len = length(hyperedge)
	num_of_edges = len * (len - 1) / 2
	I = Int[]
	J = Int[]
	V = Float64[]
	single_edge_weights = edge_weights / num_of_edges
	if len < 2
		# raise an error
		error("The number of node in an hyperedge must be at least 2.")
	else
		for i in 1:(len-1)
			for j in (i+1):len
				push!(I, hyperedge[i])
				push!(J, hyperedge[j])
				push!(V, single_edge_weights)
			end
		end
	end

	return I, J, V

end

function CliqueSampling(ar, W, mx)
	I_list = zeros(Int64, 0)
	J_list = zeros(Int64, 0)
	V_list = zeros(Float32, 0)
	for i in 1:length(ar)
		ax = length(ar[i])
		if ax <= 6
			I, J, V = CliqueExpansion(ar[i], W[i])
			append!(I_list, I)
			append!(J_list, J)
			append!(V_list, V)
		else
			# sample_num = ceil(Int32, 0.02 * ax)  # e.g., 2% of size
			# sample_num = max(sample_num, 1)      # Ensure at least 2
			I, J, V = CliqueSampler(ar[i], W[i]; sample_num = 2)
			append!(I_list, I)
			append!(J_list, J)
			append!(V_list, V)
		end
	end

	adj = sparse(I_list, J_list, V_list, mx, mx)
	AC = adj + adj'
	return AC
end


## Input: a set of random vectors, smoothing steps, star matrix, number of nodes in hypergraph
# index of the first selected smoothed vector, interval among the selected smoothed vectors, total number of smoothed vectors
# Output: a set of smoothed vectors
function Filter(rv, k, AD, mx, initial, interval, Ntot)

	sz = size(AD, 1)

	V = zeros(mx, Ntot)

	sm_vec = zeros(mx, k)

	AD = AD .* 1.0

	AD[diagind(AD, 0)] = AD[diagind(AD, 0)] .+ 0.1

	dg = sum(AD, dims = 1) .^ (-0.5)

	I2 = 1:sz

	D = sparse(I2, I2, sparsevec(dg))

	on = ones(Int, length(rv))

	sm_ot = rv - ((dot(rv, on) / dot(on, on)) * on)

	sm = sm_ot ./ norm(sm_ot)

	count = 1

	for loop in 1:k

		sm = D * sm

		sm = AD * sm

		sm = D * sm

		sm_ot = sm - ((dot(sm, on) / dot(on, on)) * on)

		sm_norm = sm_ot ./ norm(sm_ot)

		sm_vec[:, loop] = sm_norm[1:mx]

	end # for loop

	V = sm_vec[:, interval:interval:end]

	return V

end #end of function

## Input: hypergraph array, and a set of smoothed vectors
# Output: hyperedge scores
function HSC(ar, SV)
	score = zeros(eltype(SV), length(ar))
	@inbounds Threads.@threads for i in eachindex(ar)
		nodes = ar[i]
		SS = axes(SV, 2)
		for j in axes(SV, 2)
			mx, mn = -Inf, +Inf
			for node in nodes
				x = SV[node, j]
				mx = ifelse(x > mx, x, mx)
				mn = ifelse(x < mn, x, mn)
			end
			score[i] += (mx - mn)^2
		end
	end
	return score
end

# Optimized: Find neighbors only for isolated clusters
function build_isolated_cluster_neighbors(ar_new::Vector{<:AbstractVector{<:Integer}},
	isolated_cluster_ids::Vector{Int})

	# Only allocate for isolated clusters
	isolated_set = Set(isolated_cluster_ids)
	nbrs = Dict{Int, Set{Int}}()

	# Initialize empty neighbor sets for isolated clusters
	for iso_id in isolated_cluster_ids
		nbrs[iso_id] = Set{Int}()
	end

	# Process each hyperedge
	for he in ar_new
		length(he) < 2 && continue

		# Check if this hyperedge contains any isolated clusters
		isolated_in_he = [c for c in he if c in isolated_set]
		isempty(isolated_in_he) && continue

		# For each isolated cluster in this hyperedge, 
		# add all other clusters as neighbors
		for iso_c in isolated_in_he
			for other_c in he
				if iso_c != other_c
					push!(nbrs[iso_c], other_c)
				end
			end
		end
	end

	# Convert to the format expected by your merging code
	# Return a dict mapping isolated_cluster_id -> neighbors
	return Dict(iso => collect(neighbors) for (iso, neighbors) in nbrs)
end

## Input: hypergraph array
# Output: an array showing the hyperedges belong to each node
function HyperNodes(ar)

	H = INC(ar)

	NH1 = Any[]

	rr1 = H.rowval

	cc1 = H.colptr

	for i ∈ 1:size(H, 2)

		st = cc1[i]

		ed = cc1[i+1] - 1

		push!(NH1, rr1[st:ed])

	end

	return NH1

end

## The output file shows the cluster that every nodes belong to it
function CLidx(idx_mat)
	V = 1:maximum(idx_mat[end])
	for ii ∈ 1:length(idx_mat)
		idx1 = idx_mat[end-ii+1]
		idx1 = filter(index -> index <= length(V), idx1)
		if isempty(idx1)  # Check if idx1 becomes empty to avoid errors
			break
		end
		V = V[idx1]
	end
	return V
end
# end # module