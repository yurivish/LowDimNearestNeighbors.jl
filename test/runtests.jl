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

# println(SSS.satadd([0x01, 0xdd], 100))
# println(SSS.satsub([0x01, 0xdd], 100))

#

sqdist = SSS.sqdist
function nearest_linear(arr, q)
	# arr[indmin(p -> sqdist(p, q), arr)]
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
		println("Wtf. Searching for ", q, ":")
		println("\tResult: ", result, "\t", d_result)
		println("\tLinear: ", result_linear, "\t", d_linear)
		println("\t% error: ", 100 * (1 - d_linear / d_result))
		println()
	end
end

pts = [[rand(Uint8), rand(Uint8)] for i in 1:1000]
preprocess!(pts)
qs = [[rand(Uint8), rand(Uint8)] for i in 1:1000]
for q in qs
	result = nearest(pts, q)
	check(result, pts, q)
end

immutable Vec2
	x::Uint8
	y::Uint8
end
Base.getindex(v::Vec2, n::Int) = n == 1 ? v.x : n == 2 ? v.y : throw("Vec2 indexing error.")
Base.length(v::Vec2) = 2
Base.rand(::Type{Vec2}) = Vec2(rand(Uint8), rand(Uint8))

immutable Vec3
	x::Uint8
	y::Uint8
	z::Uint8
end
Base.getindex(v::Vec3, n::Int) = n == 1 ? v.x : n == 2 ? v.y : n == 3 ? v.z : throw("Vec3 indexing error.")
Base.length(v::Vec3) = 3
Base.rand(::Type{Vec3}) = Vec3(rand(Uint8), rand(Uint8), rand(Uint8))

immutable Vec4
	x::Uint8
	y::Uint8
	z::Uint8
	w::Uint8
end
Base.getindex(v::Vec4, n::Int) = n == 1 ? v.x : n == 2 ? v.y : n == 3 ? v.z : n == 4 ? v.w : throw("Vec4 indexing error.")
Base.length(v::Vec4) = 4
Base.rand(::Type{Vec4}) = Vec4(rand(Uint8), rand(Uint8), rand(Uint8), rand(Uint8))

function benchmark()
	arr = [rand(Vec3) for i in 1:100000]
	sort!(arr, lt=shuffless)
	for i in 1:10
		queries = [rand(Vec3) for i in 1:100000]
		@time for q in queries
			result = nearest(arr, q, 0.0)
			# check(result, arr, q)
		end
	end
end

benchmark()