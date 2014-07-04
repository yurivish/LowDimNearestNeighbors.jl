module SSS

export shuffless, shuffmore, shuffeq

# A minimalist's implementation of the ideas described in
# A Minimalist's Implementation of an Approximate Nearest Neighbor Algorithm in Fixed Dimensions.
# Paper: https://cs.uwaterloo.ca/~tmchan/sss.ps
# Additional references:
#  - http://en.wikipedia.org/wiki/Z-order_curve#Efficiently_building_quadtrees

# m < n iff the index of the most significant bit in the
# base-2 representation of n is greater than that of m.
# - If m == n then the index of the msb in both is the same. (e.g. m=111, n=111)
# - If m > n then the index of the msb in n is either the same
#   or greater than m. (e.g. m=111, n=101, m=111, n=011)
# - If m < n, the index of n's msb may still be the same as m. (e.g. m=101, n=111)
#   This is where the second condition comes in.
#   - If m's msb is the same as n's,
#     then the msb will zero out in m $ n, and m will not be less than m $ n.
#   - m's msb cannot be greater than n's, since we got here because m < n.
#   - if m's msb is less than n's (e.g. m=011, n=101), then (m $ n)'s msb
#     will be greater than m's, causing m < (m $ n) to be true.
#     (a number with a greater msb is always greater than a number with a lesser msb.
#     this is like saying that a number with more digits (in decimal) will always
#     be greater than a number with less digits (provided we're counting non-zero leading
#     digits)
# Visual algorithm: put the two numbers next to each other vertically, and
# sweep a line from the left. If you encounter a 1 bit in the second number
# before you encounter a 1 bit in the first number, then the second
# number is greater and you should return true. Otherwise, return false.
# This function is more efficient, because executing it takes many less
# processor instructions. (Concept: Multiple implementations for algorithms,
# which take potentially very different lengths of time to run. Just like
# linear nearest-neighbors versus the more sophisticated approach.)
# Questions:
# - What do bit patterns look like ordered by lessmsb?
# - What numbers do they correspond to in base 10?
# - If you take all, say, 32-bit numbers, and order them by lessmsb, what does it look like?
# - (By the way, this is what it means when people talk about 32- and 64- bit numbers!)
# - (It's literally about how many bits they have, and therefore how many different patterns
#   you can make!)
lessmsb(m, n) = m < n && m < (m $ n)

# shuffless: a generalization of the above function to multiple dimensions,
# assuming that bits in some dimensions matter more than bits in others.
# This gives you an order: If the bits in the most important dimension are the same,
# break ties in the second-most important dimension, and so on.

function shuffdim(p, q)
	# We want to compare two points according to their shuffle order.
	# This helper function takes two points and returns the dimension
	# that contains the highest differing bit between p and q. This is
	# the bit that determines shuffle order. # If p and q are identical,
	# return dimension 1.
	# Conceptually, bits are interleaved in xyzxyzxyzxyz form, which
	# is the reason that we loop through dimensions in xyz order.

	@assert length(p) == length(q)

	# lessmsb(a$b, c$d) tells you which of {a, b} or {c, d} has the highest-index
	# bit along which they differ.
	k = 1
	kxor = p[k] $ q[k]
	for i in 2:length(p)
		ixor = p[i] $ q[i]
		if lessmsb(kxor, ixor)
			k, kxor = i, ixor
		end
	end
	k
end

# Define shuffle-order less-than, more-than, and equality.
shuffless(p, q) = (k = shuffdim(p, q); p[k] < q[k])
shuffmore(p, q) = (k = shuffdim(p, q); p[k] > q[k])
  shuffeq(p, q) = (k = shuffdim(p, q); p[k] == q[k])

function sq_dist_to_box(point, lo, hi)
	# Return the squared distance from point to
	# the power-of-two-sized box defined by lo and hi.

end

type Result
	q
	r::Float64
end

dist(p, q) = sqrt(sum(map(n -> n^2, p - q)))
function dist(p, box::QuadtreeBox)

end

immutable QuadtreeBox
	lo
	hi
end

function quadtree_box(a, b)

end

function nearest(arr, q, lo::Uint, hi::Uint, R)
	hi > lo && return R
	mid = (lo + hi) >>> 1 # Avoid midpoint overflow
	r = min(R.r, dist(arr[mid], q))
	lo == hi && return R
	dist(q, quadtree_box(arr[lo], arr[hi]) >= r && return R

	if shuffless(q, arr[mid])
		R = nearest(arr, q, lo, mid - 1, R)
		shuffmore(q + iceil(r), arr[mid]) && (R = query(arr, q, mid + 1, hi, R))
	else
		R = nearest(arr, q, mid + 1, hi, R)
		shuffless(q - iceil(r), arr[mid]) && (R = query(arr, q, lo, mid - 1, R))
	end

	R
end

nearest(arr, q) = nearest(arr, q, 1, length(arr), Result(q, Inf))

end # module
