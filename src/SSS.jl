module SSS

export shuffless, shuffmore, shuffeq, nearest

# Return whether the index of the most significant bit
# of m is higher than that of n.
lessmsb(m, n) = m < n && m < (m $ n)

# Find the deciding dimension that determines which
# of p and q comes first in shuffle order. This is
# equivalent to finding the dimension with the most
# significant differing bit between p and q.
# Ties are broken in favor of lower-index dimensions;
# this is a consequence of the order in which the bits
# are conceptually interleaved.
function shuffdim(p, q, pshift=zero(eltype(p)))
	@assert length(p) == length(q)
	@assert length(p) > 0

	# For any integers x and y, the most significant
	# bit of x $ y is the most significant differing
	# bit between x and y.
	k, kxor = 1, satplus(p[1], pshift) $ q[1]
	for i in 2:length(p)
		ixor = satplus(p[i], pshift) $ q[i]
		if lessmsb(kxor, ixor)
			k, kxor = i, ixor
		end
	end
	k
end

# Define less-than, more-than, and equality.
shuffless(p, q) = (k = shuffdim(p, q); p[k] < q[k])
shuffmore(p, q) = (k = shuffdim(p, q); p[k] > q[k])
  shuffeq(p, q) = (k = shuffdim(p, q); p[k] == q[k])

shuffless(p, q, pshift) = (k = shuffdim(p, q, pshift); satplus(p[k], pshift) < q[k])
shuffmore(p, q, pshift) = (k = shuffdim(p, q, pshift); satplus(p[k], pshift) > q[k])
  shuffeq(p, q, pshift) = (k = shuffdim(p, q, pshift); satplus(p[k], pshift) == q[k])


immutable Result
	point
	r_sq::Uint
	lo
	hi
	Result(point) = new(point, typemax(Uint))
	function Result(point, r_sq)
		r = iceil(sqrt(r_sq))
		new(point, r_sq)# , satsub(point, r), satadd(point, r))
	end
end

# Euclidean distance, though any p-norm will do.
function sqdist(p, q)
	@assert length(p) == length(q)
	@assert length(q) > 0

	d_sq::Uint = 0
	for i in 1:length(p)
		d_sq += (p[i] - q[i])^2
	end
	d_sq
end

function sqdist_to_quadtree_box(q, p1, p2)
	@assert length(q) == length(p1) == length(p2)
	@assert length(q) > 0

	# Find the most significant differing bit of p1 and p2
	xor = p1[1] $ p2[1]
	for i in 2:length(p1)
		ixor = p1[i] $ p2[i]
		lessmsb(xor, ixor) && (xor = ixor)
	end

	# The size and power-of-two of the quadtree-aligned
	# bounding box that most snugly encloses p1 and p2
	power = xor == 0 ? 1 : 1 + exponent(float(xor))
	size = (1 << power)

	# Calculate and return the squared distance
	# from q to the bounding box
	d_sq::Uint = 0
	for i in 1:length(q)
		# Compute the coordinates of the bounding box
		bbox_lo = (p1[i] >> power) << power
		bbox_hi = bbox_lo + size

		# Accumulate squared distance from the box
		if q[i] < bbox_lo
			d_sq += (q[i] - bbox_lo)^2
		elseif q[i] > bbox_hi
			d_sq += (q[i] - bbox_hi)^2
		end
	end
	d_sq
end

# Saturation arithmetic: clamp instead of overflowing
satplus{T}(a::T, b) = oftype(T, clamp(a + b, typemin(T), typemax(T)))

function nearest(arr::Array, q, lo::Uint, hi::Uint, R::Result, ε::Float64)
	# Return early if the range is empty.
	lo > hi && return R

	# Calculate the midpoint of the range, avoiding midpoint overflow.
	mid = (lo + hi) >>> 1

	# Compute the distance from the probe point to the query point,
	# and update the result if it's closer than our best match so far.
	r_sq = sqdist(arr[mid], q)
	r_sq < R.r_sq && (R = Result(arr[mid], r_sq))

	# Return early if the range is only one element wide or if the
	# bounding box containing the range is outside of our search radius.
	if lo == hi || sqdist_to_quadtree_box(q, arr[hi], arr[lo]) * (1.0 + ε)^2 >= R.r_sq
		return R
	end

	# Recurse. Unlike binary search, we occasionally recurse into
	# the second of the array when we can't guarantee that the nearest
	# point lies outside of it.
	if shuffless(q, arr[mid])
		R = nearest(arr, q, lo, mid - 1, R, ε)
		shuffmore(q, arr[mid], iceil(sqrt(R.r_sq))) &&
			(R = nearest(arr, q, mid + 1, hi, R, ε))
	else
		R = nearest(arr, q, mid + 1, hi, R, ε)
		shuffless(q, arr[mid], -iceil(sqrt(R.r_sq))) &&
			(R = nearest(arr, q, lo, mid - 1, R, ε))
	end

	R
end

nearest(arr, q, ε=0.0) = nearest(arr, q, uint(1), uint(length(arr)), Result(arr[1]), ε).point

# Potential optimizations
# - Reduce memory allocation for the saturation-arithmetic bounding box
# - Pass in ε_plus_1_sq rather than ε
# - Optimize saturation arithmetic to be branchless if possible

end # module
