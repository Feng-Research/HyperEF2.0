using Statistics
using Base.Threads
using Base.Threads: atomic_add!

function HyperEF_2(ar, L, R, name)

	ar_new = Any[]

	idx_mat = Any[]

	Neff = zeros(Float64, mxF(ar))

	W = ones(Float64, length(ar))

	@inbounds for loop ∈ 1:L

		# println("Coarsening level: $loop")

		mx = mxF(ar)

		## expansions
		A = StarW(ar, W)

		AC = CliqueW(ar, W)
		# AC = CliqueSampling(ar, W, mx)


		## computing the smoothed vectors
		initial = 0
		SmS = 300
		interval = 20
		Nrv = 1
		Nsm = Int((SmS - initial) / interval)
		Ntot = Nrv * Nsm
		SV = zeros(Float64, mx, Ntot)
		for ii ∈ 1:Nrv
			sm = zeros(mx, Nsm)
			Random.seed!(1)
			randstring()
			rv = (rand(Float64, size(A, 1), 1) .- 0.5) .* 2
			sm = Filter(rv, SmS, A, mx, initial, interval, Nsm)
			SV[:, ((ii-1)*Nsm+1):(ii*Nsm)] = sm
		end

		SVC = zeros(Float64, mx, Ntot)
		for ii ∈ 1:Nrv
			sm = zeros(mx, Nsm)
			Random.seed!(1)
			randstring()
			rv = (rand(Float64, size(AC, 1), 1) .- 0.5) .* 2
			sm = Filter(rv, SmS, AC, mx, initial, interval, Nsm)
			SVC[:, ((ii-1)*Nsm+1):(ii*Nsm)] = sm
		end


		SVp = zeros(Float64, mx, 2 * Ntot)
		## Smoothed vectors Pool
		SVp = hcat(SV, SVC)

		## Make all the smoothed vectors orthogonal to each other
		QRpool = qr(SVp)

		SVp = Matrix(QRpool.Q)

		Q_H = zeros(Float64, size(SVp, 2))

		## Computing the ratios
		Eratio_p = zeros(Float64, length(ar), size(SVp, 2))
		for jj ∈ 1:size(SVp, 2)
			hscore = HSC(ar, SVp[:, jj])
			Q_H[jj] = sum(hscore)
			Eratio_p[:, jj] = hscore ./ sum(hscore)
		end

		# ## Approximating the effective resistance of hyperedges
		# ## by selecting the top ratio
		E2 = sort(Eratio_p, dims = 2, rev = true)
		Evec = E2[:, 1]


		# Adding the effective resistance of super nodes from previous levels
		@inbounds for kk ∈ 1:length(ar)
			nd2 = ar[kk]
			Evec[kk] = Evec[kk] + sum(Neff[nd2])
		end


		## Normalizing the ERs
		P = Evec ./ maximum(Evec)

		##
		RT = R * maximum(Evec)


		## Choosing a ratio of the hyperedges for contraction
		Nsample = length(findall(x -> x <= RT, Evec))


		PosP = sortperm(P[:, 1])
		W[PosP[1:Nsample]] = W[PosP[1:Nsample]] .* (1 .+ 1 ./ P[PosP[1:Nsample]])

		Pos = sortperm(W, rev = true)
		## low-ER diameter clustering which starts by contracting
		# the hyperedges with low ER diameter
		flag = falses(mx)

		val = 1

		idx = zeros(Int, mx)

		Neff_new = zeros(Float64, 0)

		@inbounds for ii ∈ 1:Nsample

			nd = ar[Pos[ii]]

			fg = flag[nd]

			fd1 = findall(x -> x == 0, fg)


			if length(fd1) > 1

				nd = nd[fd1]

				idx[nd] .= val

				# flag the discovered nodes
				flag[nd] .= 1

				val += 1

				## creating the super node weights
				new_val = Evec[Pos[ii]] + sum(Neff[nd])

				append!(Neff_new, new_val)

			end # endof if

		end #end of for ii

		## indexing the isolated nodes
		fdz = findall(x -> x == 0, idx)

		V = vec(val:(val+length(fdz)-1))

		idx[fdz] = V
		## Adding the weight of isolated nodes
		append!(Neff_new, Neff[fdz])

		## generating the coarse hypergraph
		ar_new = Any[]
		@inbounds for ii ∈ 1:length(ar)
			nd = ar[ii]
			nd_new = unique(idx[nd])
			push!(ar_new, sort(nd_new))
		end #end of for ii

		## Keeping the edge weights of non unique elements
		fdnu = unique(z -> ar_new[z], 1:length(ar_new))
		W2 = W[fdnu]

		## removing the repeated hyperedges
		ar_new = ar_new[fdnu]

		### removing hyperedges with cardinality of 1
		HH = INC(ar_new)
		ss = sum(HH, dims = 2)
		fd1 = findall(x -> x == 1, ss[:, 1])
		deleteat!(ar_new, fd1)
		deleteat!(W2, fd1)

		## Finding the cluster neighbors
		Neff_atomic = [Atomic{Float64}(w) for w in Neff_new]
		nclusters = maximum(idx)
		ar_int = [Int.(e) for e in ar_new]


		isolated_cluster_ids = unique([idx[v] for v in fdz])

		
		# Only build neighbors for isolated clusters
		isolated_nbrs = build_isolated_cluster_neighbors(ar_int, isolated_cluster_ids)


		dims      = size(SVp, 2)
		centroids = zeros(Float64, nclusters, dims)
		counts    = zeros(Int, nclusters)

		## calculating each cluster's effective resistance
		@inbounds for v in 1:length(idx)
			c               = idx[v]
			centroids[c, :] .+= SVp[v, :]
			counts[c]       += 1
		end
		centroids ./= counts

		mxd = length(idx)
		inc_weight = zeros(Float64, mxd)

		for (he_idx, he) in enumerate(ar)
			w = Evec[he_idx]
			@inbounds for v in he
				inc_weight[v] += w
			end
		end


		invQH = 1.0 ./ Q_H

		## Resistance based Local clustering
		@threads for t in eachindex(fdz)
			v = fdz[t]
			iso_id = idx[v]
			neighIDs = get(isolated_nbrs, iso_id, Int[])

			isempty(neighIDs) && continue

			rowvec = reshape(@view(SVp[v, :]), 1, :)
			Delta = abs.(centroids[neighIDs, :] .- rowvec)
			dists = vec(maximum(Delta .* reshape(invQH, 1, :), dims = 2))
			best = neighIDs[argmin(dists)]			#Find minimum distance

			idx[v] = best							#Merge with best cluster

			inc_vertex = Neff_atomic[iso_id][]
			# edge_sum = inc_weight[v]
			atomic_add!(Neff_atomic[best], inc_vertex)
		end

		#calculate new supernode effective Resistance
		Neff_new = Float64[x[] for x in Neff_atomic]


		ar_new1 = Any[]
		@inbounds for ii ∈ 1:length(ar)

			nd = ar[ii]
			nd_new = unique(idx[nd])
			push!(ar_new1, sort(nd_new))
		end 

		## Keeping the edge weights of non unique elements
		fdnu = unique(z -> ar_new1[z], 1:length(ar_new1))

		W3 = W[fdnu]

		## removing the repeated hyperedges
		ar_new1 = ar_new1[fdnu]

		### removing hyperedges with cardinality of 1
		HH = INC(ar_new1)
		ss = sum(HH, dims = 2)
		fd1 = findall(x -> x == 1, ss[:, 1])
		deleteat!(ar_new1, fd1)
		deleteat!(W3, fd1)

		Neff = Neff_new

		W = W3

		push!(idx_mat, idx)
		global idxx = CLidx(idx_mat)
		ar = ar_new1
	end #end for loop

	## output cluster file
	filename = "$name.idx"
	open(filename, "w") do io
		for ii in 1:length(idxx)
			println(io, idxx[ii])
		end
	end
end
