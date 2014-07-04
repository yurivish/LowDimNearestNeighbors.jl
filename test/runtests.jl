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

satadd = SSS.satadd
println(satadd([0xdd, 0x22], 100))

satsub = SSS.satsub
println(satsub([0xdd, 0x22], 100))

#

quadtree_box = SSS.quadtree_box
# println(quadtree_box(1, 1))
# println(quadtree_box(1, 2))
# println(quadtree_box(1, 3))
# println(quadtree_box(1, 4))
# println(quadtree_box(7, 8))
# println(quadtree_box(7, 7))
# println(quadtree_box(8, 8))
# println(quadtree_box(9, 9))
# println(quadtree_box(10, 10))
# println(quadtree_box(11, 11))
# println(quadtree_box([7, 16], [11, 18]))
for i in 1:1000
	# a = Uint8[25,59]
	# b = Uint8[37,6]
	a = [rand(Uint8), rand(Uint8)]
	b = SSS.satadd(a, uint8(10))
	@test shuffless(a, b)
	box = quadtree_box(a, b)
	println("Box for $a, $b: ", box)
	@test shuffless(box.lo, a)
	@test shuffless(a, b)
	@test shuffless(b, box.hi) || shuffeq(b, box.hi)
end



#

# dist = SSS.dist
# function nearest_linear(arr, q)
# 	local best
# 	best_dist = Inf
# 	for point in arr
# 		if dist(point, q) < best_dist
# 			best, best_dist = point, dist(point, q)
# 		end
# 	end
# 	best
# end

# pts = [[rand(Uint8), rand(Uint8)] for i in 1:1000]
# sort!(pts, lt=shuffless)
# for i in 1:10
# 	pt = [rand(Uint8), rand(Uint8)]
# 	result = nearest(pts, pt)
# 	result_linear = nearest_linear(pts, pt)
# 	if dist(pt, result) != dist(pt, result_linear)
# 		println("Nearest point to ", pt, ": ", result, "; linear=", result_linear)
# 		println("--- Distances: ", dist(pt, result), ", ", dist(pt, result_linear))
# 		println()
# 	end
# end