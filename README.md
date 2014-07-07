# Nearest-Neighbor Search in Low Dimensions

This package implements approximate nearest-neighbor search in low dimensions for points with integer coordinates, using an elegant idea from a 2006 [paper](http://cs.uwaterloo.ca/~tmchan/sss.ps) by Timothy M. Chan.

```julia
	using SSS

	# Create an array of random 3d points
	arr = [[rand(Uint8), rand(Uint8), rand(Uint8)] for i in 1:100000]

	# Preprocess it to prepare for efficient searching
	preprocess!(arr)

	# Perform an exact nearest-neighbor search
	query  = [rand(Uint8), rand(Uint8), rand(Uint8)]
	result = nearest(arr, query)
	println("Nearest point to $query: $result")
	println("Distance: ", sqrt(SSS.sqdist(query, result)))

	# Perform an approximate nearest-neighbor search to find
	# a point whose distance to the query is within 25% of 
	# the best possible result.
	result = nearest(arr, query, 0.25)
	println("Approximate nearest point to $query: $result")
	println("Distance: ", sqrt(SSS.sqdist(query, result)))
```

## Notes

The algorithm is _in-place_, i.e. it requires no extra space beyond the input array. Instead, spatial information is encoded in the permutation of points-- the preprocessing step sorts the array to prepare for efficient queries.

The approach here works best in low dimensions (2, 3, 4), but the code is generic and will work for points of arbitrary dimension so long as they implement `getindex` and `length`.

When performing approximate searches, the points found by the algorithm tend to be better than you'd expect based on the approximation factor. For example, in the program above it will often find an exact match.

Currently, there is risk of overflow in the distance function -- if the squared euclidean distance exceeds `typemax(Uint)`, overflow will occur and the results returned by the algorithm will be incorrect.

This code has been used to some success to implement nearest-neighbor search in RGB colorspace. If you use it for anything interesting, [let me know](mailto:yurivish@gmail.com)!

[![Build Status](https://travis-ci.org/yurivish/SSS.jl.svg?branch=master)](https://travis-ci.org/yurivish/SSS.jl)
