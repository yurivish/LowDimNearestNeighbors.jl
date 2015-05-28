using LowDimNearestNeighbors
using Base.Test

# Define multidimensional vector types for testing
typealias Vec2{T <: Unsigned} Tuple{T, T}
Base.call{T}(::Type{Vec2{T}}, x::T, y::T) = (unsigned(x), unsigned(y))
Base.rand{T}(::Type{Vec2{T}}) = Vec2(rand(T), rand(T))

#

lessmsb = LowDimNearestNeighbors.lessmsb
@test  lessmsb(1, 2) # msb(001) <  msb(010) 
@test !lessmsb(2, 1) # msb(010) >  msb(001) 

@test !lessmsb(2, 3) # msb(010) == msb(011)
@test !lessmsb(3, 2) # msb(011) == msb(010)

@test  lessmsb(1, 5) # msb(001) <  msb(101)
@test !lessmsb(5, 1) # msb(101) >  msb(001)

#

shuffdim = LowDimNearestNeighbors.shuffdim
# The shuffdim should be 1 for identical points
@test shuffdim([1, 2, 3], [1, 2, 3]) == 1
@test shuffdim([123, 123], [123, 123]) == 1

# The first dimension is greater in p
@test shuffdim([1, 0], [0, 0]) == 1
@test shuffdim([1, 0], [0, 1]) == 1

# The first dimension is greater in q
@test shuffdim([0, 0], [1, 0]) == 1
@test shuffdim([0, 1], [1, 0]) == 1

# The second dimension is greater in p
@test shuffdim([1, 2], [1, 1]) == 2
# The second dimension is greater in q
@test shuffdim([1, 1], [1, 2]) == 2

#

@test  shuffless([1, 1, 3], [1, 2, 3])
@test !shuffless([1, 2, 3], [1, 2, 3])
@test !shuffless([1, 3, 3], [1, 2, 3])

@test !shuffmore([1, 1, 3], [1, 2, 3])
@test !shuffmore([1, 2, 3], [1, 2, 3])
@test  shuffmore([1, 3, 3], [1, 2, 3])

@test !shuffeq([1, 3, 3], [1, 2, 3])
@test shuffeq([1, 2, 3], [1, 2, 3])

# Test preprocess!
let
	arr = [Vec2(1, 1), Vec2(0, 1), Vec2(1, 0), Vec2(0, 0)]
	preprocess!(arr)
	@test arr[1] == Vec2(0, 0)
	@test arr[2] == Vec2(0, 1)
	@test arr[3] == Vec2(1, 0)
	@test arr[4] == Vec2(1, 1)

	# Error on negative coordinates
	@test_throws ErrorException preprocess!([1, 2, -3])

	# Error on noninteger coordinates
	@test_throws ErrorException preprocess!([1, 2.2, 3])
	@test_throws ErrorException preprocess!([1, 2.0, 3])
end

# Test Shifted indexing, clamping, and length
let
	el = Vec2{UInt64}(1, 2)
	shifted = LowDimNearestNeighbors.ShiftedPos{Vec2{UInt64}}(el, 5)
	@test shifted[1] == 6
	@test shifted[2] == 7
	@test length(shifted) == 2

	el = Vec2{UInt64}(4, 6)
	shifted = LowDimNearestNeighbors.ShiftedNeg{Vec2{UInt64}}(el, 5)
	@test shifted[1] == 0
	@test shifted[2] == 1
	@test length(shifted) == 2
end

# Test sqdist
let
	sqdist = LowDimNearestNeighbors.sqdist
	@test sqdist(2, 2) == 0
	@test sqdist(2, 3) == 1
	@test sqdist(2, 4) == 4
	@test sqdist(Vec2(3, 0), Vec2(0, 4)) == 5*5

	# Test saturation cases
	@test sqdist(0, typemax(UInt)) == typemax(UInt)
	@test sqdist(typemax(UInt), 0) == typemax(UInt)
	@test sqdist(Vec2(3, 3), Vec2(typemax(UInt), typemax(UInt))) == typemax(UInt)
end

# Test sqdist_to_quadtree_box
let
	sqdist_to_quadtree_box = LowDimNearestNeighbors.sqdist_to_quadtree_box
	@test sqdist_to_quadtree_box(Vec2(0, 0), Vec2(3, 0), Vec2(0, 3)) == 0
	@test sqdist_to_quadtree_box(Vec2(2, 2), Vec2(3, 0), Vec2(0, 3)) == 0

	@test sqdist_to_quadtree_box(Vec2(0, 0), Vec2(0, 0), Vec2(3, 3)) == 0
	@test sqdist_to_quadtree_box(Vec2(2, 2), Vec2(0, 0), Vec2(3, 3)) == 0

	@test sqdist_to_quadtree_box(Vec2(0, 0), Vec2(2, 2), Vec2(3, 3)) == 8
	@test sqdist_to_quadtree_box(Vec2(4, 4), Vec2(0, 0), Vec2(1, 1)) == 8
end


# Nearest: Deterministic test
let
	arr = preprocess!([Vec2(3, 3), Vec2(10, 0), Vec2(5, 8)])
	@test nearest(arr, Vec2(0, 0)) == Vec2(3, 3)
	@test nearest(arr, Vec2(8, 1)) == Vec2(10, 0)
	@test nearest(arr, Vec2(1000, 1000)) == Vec2(5, 8)
	@test nearest(arr, Vec2(5, 8)) == Vec2(5, 8)
end

# Nearest: Randomized test
let
	sqdist = LowDimNearestNeighbors.sqdist
	function linear_nearest(arr, q)
		local best
		best_sqdist = Inf
		for point in arr
			if sqdist(point, q) < best_sqdist
				best, best_sqdist = point, sqdist(point, q)
			end
		end
		best
	end

	function test(arr, queries; verbose=false)
		for q in queries
			result = nearest(arr, q)
			result_sqdist = sqdist(q, result)

			correct_result = linear_nearest(arr, q)
			correct_sqdist = sqdist(q, correct_result)

			if verbose && result_sqdist != correct_sqdist
				result_dist = sqrt(result_sqdist)
				correct_dist = sqrt(correct_sqdist)
				println("Mismatch when searching for ", q, ":")
				println("\t Result: ", result, "\t", result_dist)
				println("\tCorrect: ", correct_result, "\t", correct_dist)
				println("\t% error: ", 100 * (1 - correct_dist / result_dist), "%")
				println()
			end

			@test result_sqdist == correct_sqdist
		end
	end

	srand(0)
	numelements = numqueries = 1000

	test(
		preprocess!([rand(Vec2{UInt64}) for i in 1:numelements]),
		[rand(Vec2{UInt64}) for i in 1:numqueries]
	)

	test(
		preprocess!([rand(Vec2{UInt32}) for i in 1:numelements]),
		[rand(Vec2{UInt32}) for i in 1:numqueries]
	)

	test(
		preprocess!([rand(Vec2{UInt16}) for i in 1:numelements]),
		[rand(Vec2{UInt16}) for i in 1:numqueries]
	)

	test(
		preprocess!([rand(Vec2{UInt8}) for i in 1:numelements]),
		[rand(Vec2{UInt8}) for i in 1:numqueries]
	)

end

# Benchmarks
let
	function benchmark(numelements, numqueries)
		for i in 1:10
			arr = preprocess!([rand(Vec2{UInt8}) for i in 1:numelements])
			queries = [rand(Vec2{UInt8}) for i in 1:numqueries]
			@time for q in queries
				nearest(arr, q)
			end
		end
	end

	# benchmark(100000, 100000)
end
