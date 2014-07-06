using SSS
using Base.Test

lessmsb = SSS.lessmsb
@test  lessmsb(1, 2) # msb(001) <  msb(010) 
@test !lessmsb(2, 1) # msb(010) >  msb(001) 

@test !lessmsb(2, 3) # msb(010) == msb(011)
@test !lessmsb(3, 2) # msb(011) == msb(010)

@test  lessmsb(1, 5) # msb(001) <  msb(101)
@test !lessmsb(5, 1) # msb(101) >  msb(001)

#

shuffdim = SSS.shuffdim
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

#

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

# Test preprocess!
let
	T = Vec2{Int64}
	arr = [T(1, 1), T(0, 1), T(1, 0), T(0, 0)]

	preprocess!(arr)
	@test arr[1] == T(0, 0)
	@test arr[2] == T(0, 1)
	@test arr[3] == T(1, 0)
	@test arr[4] == T(1, 1)
end

# Test Shifted indexing and length
let
	el = Vec2{Int64}(1, 2)
	shifted = SSS.Shifted{Vec2{Int64}}(el, 5)
	@test shifted[1] == 6
	@test shifted[2] == 7
	@test length(shifted) == 2
end

# Test sqdist
let
	sqdist = SSS.sqdist
	@test sqdist(2, 2) == 0
	@test sqdist(2, 3) == 1
	@test sqdist(2, 4) == 4
	@test sqdist(Vec2(3, 0), Vec2(0, 4)) == 5*5
end

# Test sqdist_to_quadtree_box
let
	sqdist_to_quadtree_box = SSS.sqdist_to_quadtree_box
	@test sqdist_to_quadtree_box(Vec2(0, 0), Vec2(3, 0), Vec2(0, 3)) == 0
	@test sqdist_to_quadtree_box(Vec2(2, 2), Vec2(3, 0), Vec2(0, 3)) == 0

	@test sqdist_to_quadtree_box(Vec2(0, 0), Vec2(0, 0), Vec2(3, 3)) == 0
	@test sqdist_to_quadtree_box(Vec2(2, 2), Vec2(0, 0), Vec2(3, 3)) == 0

	@test sqdist_to_quadtree_box(Vec2(0, 0), Vec2(2, 2), Vec2(3, 3)) == 8
	@test sqdist_to_quadtree_box(Vec2(4, 4), Vec2(0, 0), Vec2(1, 1)) == 8
end

# Test nearest
let
	sqdist = SSS.sqdist

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
		arr = preprocess!([rand(Vec2{Uint8}) for i in 1:numelements])
		queries = [rand(Vec2{Uint8}) for i in 1:numqueries]
		for q in queries
			result = nearest(arr, q)
			correct_result = linear_nearest(arr, q)

			distsq = sqdist(q, result)
			correct_distsq = sqdist(q, correct_result)

			if verbose && distsq != correct_distsq
				result_dist = sqrt(distsq)
				correct_dist = sqrt(correct_distsq)
				println("Mismatch when searching for ", q, ":")
				println("\t Result: ", result, "\t", result_dist)
				println("\tCorrect: ", correct_result, "\t", correct_dist)
				println("\t% error: ", 100 * (1 - correct_dist / result_dist), "%")
				println()
			end

			@test distsq == correct_distsq
		end
	end

	srand(0)
	test(1000, 10000)
end
