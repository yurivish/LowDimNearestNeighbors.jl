# Shift-Shuffle-Sort: Approximate Nearest-Neighbor Search in Low Dimensions

This package implements approximate nearest-neighbor search in multiple dimensions for points with integer coordinates. The idea is from a 2006 [paper](http://cs.uwaterloo.ca/~tmchan/sss.ps) by Timothy M. Chan.

```julia
	using SSS

	# Create an array of random points in 3d
	arr = [[rand(Uint8), rand(Uint8), rand(Uint8)] for i in 1:10000]

	# Preprocess the array to prepare it for efficient searches
	preprocess!(arr)

	query  = [100, 100]
	result = nearest(arr, query)
	println("Nearest point to $query: $result")
```

## Notes

The algorithm is _in-place_, i.e. it requires no extra space beyond the input array. Instead, spatial information is encoded in the permutation of points: the preprocessing step sorts the array in a particular way to prepare for efficient queries.

[![Build Status](https://travis-ci.org/yurivish/SSS.jl.svg?branch=master)](https://travis-ci.org/yurivish/SSS.jl)
