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

sqdist = SSS.sqdist
function nearest_linear(arr, q)
	best, best_sqdist = q, Inf
	for point in arr
		if sqdist(point, q) < best_sqdist
			best, best_sqdist = point, sqdist(point, q)
		end
	end
	best
end

function check(result, arr, q)
	result_linear = nearest_linear(arr, q)
	sqd_result = sqdist(q, result)
	sqd_linear = sqdist(q, result_linear)
	if sqd_result != sqd_linear
		d_result = sqrt(sqd_result)
		d_linear = sqrt(sqd_linear)
		println("Mismatch when searching for ", q, ":")
		println("\tResult: ", result, "\t", d_result)
		println("\tLinear: ", result_linear, "\t", d_linear)
		println("\t% error: ", 100 * (1 - d_linear / d_result))
		println()
	end
end

points = preprocess!([[rand(Uint8), rand(Uint8)] for i in 1:1000])
qs = [[rand(Uint8), rand(Uint8)] for i in 1:1000]
for q in qs
	result = nearest(points, q)
	check(result, points, q)
end

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

function benchmark()
	arr = preprocess!([rand(Vec3{Uint8}) for i in 1:100000])
	for i in 1:10
		queries = [rand(Vec3{Uint8}) for i in 1:100000]
		@time for q in queries
			result = nearest(arr, q, 0.0)
			# check(result, arr, q)
		end
	end
end

benchmark()