# Nearest-Neighbor Search in Low Dimensions

This package implements approximate nearest-neighbor search in low dimensions for points with integer coordinates, using an elegant idea from a 2006 [paper](http://cs.uwaterloo.ca/~tmchan/sss.ps) by Timothy M. Chan.

```julia
	using SSS

	# Create an array of random points in 3d
	arr = [[rand(Uint8), rand(Uint8), rand(Uint8)] for i in 1:10000]

	# Preprocess the array to prepare it for efficient searching
	preprocess!(arr)

	# Perform an exact nearest-neighbor search
	query  = [100, 100, 100]
	result = nearest(arr, query)
	println("Nearest point to $query: $result")

	# Perform an approximate nearest-neighbor search
	# to find a point whose distance to the query is
	# within 15% of the best possible result.
	result = nearest(arr, query, 0.15)
	println("Approximate nearest point to $query: $result")
```

## Notes

The approach here works best in low dimensions (2, 3, 4), but the code is generic and will work for points of arbitrary dimension so long as they implement `getindex` and `length`.

The algorithm is _in-place_, i.e. it requires no extra space beyond the input array. Instead, spatial information is encoded in the permutation of points. The preprocessing step sorts the array to prepare for efficient queries.

When performing approximate searches, the points found by the algorithm tend to be better than you'd expect based on the approximation factor.

Currently, there is risk of overflow in the distance function -- if the squared euclidean distance exceeds `typemax(Uint)`, overflow will occur and the results returned by the algorithm will be incorrect.

This code has been used to some success to implement nearest-neighbor search in RGB colorspace. If you use it for anything interesting, [let me know](mailto:yurivish@gmail.com)!

[![Build Status](https://travis-ci.org/yurivish/SSS.jl.svg?branch=master)](https://travis-ci.org/yurivish/SSS.jl)
