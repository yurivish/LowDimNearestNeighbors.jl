# Nearest-Neighbor Search in Low Dimensions

This package implements approximate nearest-neighbor search in multiple dimensions for points with integer coordinates. The idea is from a 2006 [paper](http://cs.uwaterloo.ca/~tmchan/sss.ps) by Timothy M. Chan.

```julia
	using SSS

	# Create an array of random points in 3d
	arr = [[rand(Uint8), rand(Uint8), rand(Uint8)] for i in 1:10000]

	# Preprocess the array to prepare it for efficient searches
	preprocess!(arr)

	query  = [100, 100]

	# Perform an exact nearest-neighbor search
	result = nearest(arr, query)
	println("Nearest point to $query: $result")

	# Perform an approximate nearest-neighbor search
	# to find a point whose distance to the query is
	# within 15% of the exact match
	result = nearest(arr, query, 0.15)
	println("Approximately nearest point to $query: $result")
```

## Notes

The algorithm is _in-place_, i.e. it requires no extra space beyond the input array. Instead, spatial information is encoded in the permutation of points. The preprocessing step sorts the array to prepare for efficient queries.

When performing approximate searches, the points found by the algorithm tend to be better than you'd expect based on the approximation factor.

The algorithm guarantees O(_n_ log(_n_)) preprocessing time and O(log(_n_)) query time, but the primary goal of this approach is simplicity rather than optimal speed in practice.

Currently, there is risk of overflow in the distance function -- if the squared euclidean distance exceeds `typemax(Uint)`, overflow will occur and the results returned by the algorithm will be incorrect.

[![Build Status](https://travis-ci.org/yurivish/SSS.jl.svg?branch=master)](https://travis-ci.org/yurivish/SSS.jl)
