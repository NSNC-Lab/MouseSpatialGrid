using Optim, Distributed, MAT, Distributions, DelimitedFiles
# using Distances
# include("distMatrices.jl")
# function density_est(labels::Array{T,1} where {T<:Any}, data::Array{T,1} where {T<:Any}; 
#     dist_metric=Euclidean(), all = false)
#     dists = distMat(dist_metric, data)
#     return density_est(labels, dists, all=all)
# end

# Returns the bias in estimate in h points from N, for classSizes
function hypergeometricBias(N, h, classSizes)
    I = 0
    n_s = length(classSizes)
    for n_c in classSizes
        dst = Hypergeometric(n_c-1, N-n_c, h-1)
        for r=1:h
            I += pdf(dst, r-1) * log2(n_s*r/h) * (n_c/N)
        end
    end
    return I
end

# Return the estimate of MI for labels labels and distance matrix dists
# Labels : a 1d array of any type
#  Dists : a 2d array of numbers
function density_est(labels::Array{T, 1} where {T<:Any}, dists::Array{T, 2} where {T<:Number};
        all = false)
    C = Dict{eltype(labels), Array{Int64, 1}}()
    for c in unique(labels)
        C[c] = findall((a->a==c), labels)
    end
    CLens = map(length, values(C))
    N = length(labels)
    n_s = length(keys(C))

    upper_limit = N
    if upper_limit <= 3 return 0 end

    sorted = zeros(size(dists))
    for i in 1:N
        sorted[i, 1:Int(upper_limit)] = partialsort(view(dists, i, :), 1:Int(upper_limit))
    end

    function estimator(h)
        h = Int64(round(h))
        MI = 0
        for i in 1:N
            c = labels[i]
            bRadius = sorted[i, h]
            bContents = count((a->a<bRadius), view(dists, i, C[c]))
            bContents += (h - count((a->a<bRadius), view(dists, i, :))) * count((a->a==bRadius), view(dists, i, C[c])) / count((a->a==bRadius), view(dists, i, :))
            MI += log2(n_s * bContents / h) / N
        end
        if all return MI, hypergeometricBias(N, h, CLens) end
        return - (MI - hypergeometricBias(N, h, CLens))
    end

    if all 
        res = pmap(estimator, 1:N)
        mis = [x[1] for x in res]
        biases = [x[2] for x in res]
        return mis, biases
    end
    res = optimize(estimator, 2, N, abs_tol=1, method=GoldenSection())
    return - res.minimum
end


if length(ARGS) == 1



	#Just for testing, going to have this run 20 times to see if it is faster just to handle things within Julia
	#Looks like this might natively take a matrix. Just try redoing the distance calculations and feeding in the matrix. If this does not work then 	#you can just do the zz for loop
    #for zz = 1:20
	
    	fileLoc = ARGS[1]
    	if fileLoc[end-2:end] == "mat"

        	println("Loading $fileLoc")
        	vars = matread(fileLoc)

		#println(vars)

        	allLabels = vars["labels"]
		
		#Changing to distMat for now because that seems to be the new name for whatever reason
			
        	allDists = vars["distMat"]
        	println("Found $(size(allDists)[1]) distances, size $(size(allDists)[2:3])")
        	println("Found $(size(allLabels)[1]) labels, size $(size(allLabels)[2])")

		println(allDists)

        	n = min(size(allDists)[1], size(allLabels)[1])
        	MIs = []
		
		println("Marking here to make sure I can find n")
		println(n)
		println("spacing")

        	#for i=1:n
            	#	MI = density_est(allLabels[i,:], allDists[i,:,:])
            	#	println("$i $MI")
            	#	push!(MIs, MI)
        	#end

		#Changing to handle matrix
		for i=1:n
            	 	MI = density_est(allLabels[i,:], allDists[i,:,:])
            		println("$i $MI")
            		push!(MIs, MI)
        	end

		println("Made it here")


        	writedlm(string(chop(fileLoc,head=0,tail=14),"_MIs.csv"), MIs, ',')
    	else
        	println("Not a '.mat' file")
    	end
 
     #end

elseif length(ARGS) == 0
    println("No file location passed")
end