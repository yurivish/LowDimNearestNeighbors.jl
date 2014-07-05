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
function shuffdim(p, q)
	@assert length(p) > 0 && length(q) > 0
	@assert length(p) == length(q)

	# For any integers x and y, the most significant
	# bit of x $ y is the most significant differing
	# bit between x and y.
	k, kxor = 1, p[1] $ q[1]
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

type Result
	point
	r_sq::Float64
end

# Euclidean distance, though any p-norm will do.
# dist(p, q) = norm(int(p) - int(q))
# Squared euclidean distance, though any p-norm will do.
sqdist(p, q) = sum(map(n -> n^2, int(p) - int(q)))

function sqdist_to_quadtree_box(q, p1, p2)
	@assert length(q) == length(p1) == length(p2)
	@assert length(q) > 0 && length(p1) > 0 && length(p2) > 0

	# Find the most significant differing bit between p1 and p2
	xor = p1[1] $ p2[1]
	for i in 2:length(p1)
		ixor = p1[i] $ p2[i]
		lessmsb(xor, ixor) && (xor = ixor)
	end

	# The size and power-of-two of the smallest quadtree-aligned
	# bounding box that encompasses p1 and p2
	power = xor == 0 ? 1 : 1 + exponent(float(xor))
	size = (1 << power)

	# Calculate the distance from q to the bounding box
	sqdist = 0
	for i in 1:length(q)
		# Compute the coordinates of the bounding box in this dimension
		bbox_lo = (p1[i] >> power) << power
		bbox_hi = bbox_lo + size

		# Accumulate squared distance in this dimension
		if q[i] < bbox_lo
			sqdist += (q[i] - bbox_lo)^2
		elseif q[i] > bbox_hi
			sqdist += (q[i] - bbox_hi)^2
		end
	end
	sqdist
end

# Saturation arithmetic for shifting points without running into overflows
satadd(p, r) = map(c -> oftype(c, min(c + r, typemax(c))), p)
satsub(p, r) = map(c -> oftype(c, max(c - r, typemin(c))), p)

function nearest(arr, q, lo::Uint, hi::Uint, R, ε::Float64)
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
	(lo == hi || sqdist_to_quadtree_box(q, arr[hi], arr[lo]) * (1.0 + ε)^2 >= R.r_sq) && return R

	# Recurse, à la binary search. Unlike binary search, we occasionally recurse
	# into the second of the array when we can't guarantee that the nearest point
	# lies inside it.
	if shuffless(q, arr[mid])
		R = nearest(arr, q, lo, mid - 1, R, ε)
		shuffmore(satadd(q, iceil(sqrt(R.r_sq))), arr[mid]) && (R = nearest(arr, q, mid + 1, hi, R, ε))
	else
		R = nearest(arr, q, mid + 1, hi, R, ε)
		shuffless(satsub(q, iceil(sqrt(R.r_sq))), arr[mid]) && (R = nearest(arr, q, lo, mid - 1, R, ε))
	end

	R
end

nearest(arr, q) = nearest(arr, q, uint(1), uint(length(arr)), Result(q, Inf), 0.0).point
nearest(arr, q, ε) = nearest(arr, q, uint(1), uint(length(arr)), Result(q, Inf), ε).point

end # module
