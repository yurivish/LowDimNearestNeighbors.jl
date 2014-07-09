using LowDimNearestNeighbors
using Base.Test

# Define multidimensional vector types for testing
immutable Vec2{T}
	x::T
	y::T
end
Base.getindex(v::Vec2, n::Int) = n == 1 ? v.x : n == 2 ? v.y : throw("Vec2 indexing error.")
Base.length(v::Vec2) = 2
Base.rand{T}(::Type{Vec2{T}}) = Vec2(rand(T), rand(T))

immutable Vec3{T}
	x::T
	y::T
	z::T
end
Base.getindex(v::Vec3, n::Int) = n == 1 ? v.x : n == 2 ? v.y : n == 3 ? v.z : throw("Vec3 indexing error.")
Base.length(v::Vec3) = 3
Base.rand{T}(::Type{Vec3{T}}) = Vec3(rand(T), rand(T), rand(T))

immutable Vec4{T}
	x::T
	y::T
	z::T
	w::T
end
Base.getindex(v::Vec4, n::Int) = n == 1 ? v.x : n == 2 ? v.y : n == 3 ? v.z : n == 4 ? v.w : throw("Vec4 indexing error.")
Base.length(v::Vec4) = 4
Base.rand{T}(::Type{Vec4{T}}) = Vec4(rand(T), rand(T), rand(T), rand(T))

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

	# # Error on negative coordinates
	@test_throws ErrorException preprocess!([1, 2, -3])

	# # Error on noninteger coordinates
	@test_throws ErrorException preprocess!([1, 2.2, 3])
	@test_throws ErrorException preprocess!([1, 2.0, 3])
end

# Test Shifted indexing and length
let
	el = Vec2{Int64}(1, 2)
	shifted = LowDimNearestNeighbors.Shifted{Vec2{Int64}}(el, 5)
	@test shifted[1] == 6
	@test shifted[2] == 7
	@test length(shifted) == 2
end

# Test sqdist
let
	sqdist = LowDimNearestNeighbors.sqdist
	@test sqdist(2, 2) == 0
	@test sqdist(2, 3) == 1
	@test sqdist(2, 4) == 4
	@test sqdist(Vec2(3, 0), Vec2(0, 4)) == 5*5

	# This does not yet work:
	println(sqdist(0, typemax(Uint)))
	@test sqdist(0, typemax(Uint)) == typemax(Uint)
	@test sqdist(Vec2(3, 3), Vec2(typemax(Uint), typemax(Uint))) == typemax(Uint)
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
	@test nearest(arr, Vec2(5.7, 8.1)) == Vec2(5, 8)
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

	function test(numelements, numqueries; verbose=false)
		arr = preprocess!([rand(Vec2{Uint16}) for i in 1:numelements])
		queries = [rand(Vec2{Uint16}) for i in 1:numqueries]
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
	test(1000, 1000, verbose=false)
end

# Benchmarks
let
	function benchmark(numelements, numqueries)
		for i in 1:10
			arr = preprocess!([rand(Vec2{Uint8}) for i in 1:numelements])
			queries = [rand(Vec2{Uint8}) for i in 1:numqueries]
			@time for q in queries
				nearest(arr, q)
			end
		end
	end

	benchmark(100000, 100000)
end
