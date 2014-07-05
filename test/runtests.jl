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

pts = [[rand(Uint8), rand(Uint8)] for i in 1:1000]
sort!(pts, lt=shuffless)
qs = [[rand(Uint8), rand(Uint8)] for i in 1:1000]
for q in qs
	result = nearest(pts, q)
	result_linear = nearest_linear(pts, q)

	if sqdist(q, result) != sqdist(q, result_linear)
		println("Wtf. Searching for ", q, ":")
		println("\tResult: ", result, "\t", sqdist(result, q))
		println("\tLinear: ", result_linear, "\t", sqdist(result_linear, q))
		println()
	end
end

for i in 1:10
	qs = [[rand(Uint8), rand(Uint8)] for i in 1:10000]
	@time for q in qs
		result = nearest(pts, q)
	end
end