using Statistics
using SparseArrays
using Random
using LinearAlgebra


include("HyperEF_2_0.jl")
include("Functions.jl")



filename = "ibm01.hgr"

cd("../data/")  # Ensure this path is correct relative to Run_HyperEF.jl's location
ar = ReadInp(filename)
cd("../src/")  # Ensure this path is correct

# ----------------------------
# Set Parameters
# ----------------------------

## L: the number of coarsening levels, e.g., 1, 2, 3, 4, ...
L = 3

## R: Effective resistance threshold for growing the clusters (0 < R <= 1)
R = 1.0  # Use Float64 for precision


HyperEF_2(ar, L, R, filename)
