# An implementation of the ideas described in A Minimalist's
# Implementation of an Approximate Nearest Neighbor Algorithm
# in Fixed Dimensions.
# Paper: http://cs.uwaterloo.ca/~tmchan/sss.ps
# Further reading: http://en.wikipedia.org/wiki/Z-order_curve

module LowDimNearestNeighbors

export shuffless, shuffmore, shuffeq, preprocess!, nearest

# Return whether the index of the most significant bit
# of m is higher than that of n.
lessmsb(m, n) = m < n && m < (m $ n)

# Find the deciding dimension that determines which
# of p and q comes first in shuffle order. This is
# equivalent to finding the dimension with the most
# significant differing bit between p and q.
# Ties are broken in favor of lower-index dimensions.
function shuffdim(p, q)
	@assert length(p) == length(q)
	@assert length(p) > 0

	# For any integers x and y, the msb of x $ y is
	# the most significant differing bit between
	# x and y.
	k, kxor = 1, p[1] $ q[1]
	for i in 2:length(p)
		ixor = p[i] $ q[i]
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

# Sort the given array in shuffle order to prepare it for nearest-neighbor queries.
function preprocess!(arr)
	for p in arr
		for i in length(p)
			!(typeof(p[i]) <: Unsigned) && throw(ErrorException("All coordinates must be unsigned integers."))
		end
	end
	sort!(arr, lt=shuffless)
end

# Code for branch-free saturation arithmetic from
# http://locklessinc.com/articles/sat_arithmetic/
satadd{T <: Unsigned}(x::T, y::T) = (res = x + y; res | unsigned(-signed(T(res < x))))
satsub{T <: Unsigned}(x::T, y::T) = (res = x - y; res & unsigned(-signed(T(res <= x))))

function satmul(x::UInt32, y::UInt32)
	res = UInt64(x) * UInt64(y)
	hi::UInt32 = res >> 32
	lo::UInt32 = res & 0x0000ffff
	lo | unsigned(-Int32(hi != 0))
end

function satmul(x::UInt64, y::UInt64)
	res = UInt128(x) * UInt128(y)
	hi::UInt64 = res >> 64
	lo::UInt64 = res & 0x0000000000000000ffffffffffffffff
	lo | unsigned(-Int64(hi != 0))
end

# There's no way for the result to underflow or overflow.
satdiv{T <: Unsigned}(x::T, y::T) = x / y

# Represent shifted points by their own type.
# It is assumed that data is a point with
# eltype(data) <: Unsigned and the value
# of offset fits into the type eltype(data).
immutable ShiftedPos{Q}
	data::Q
	offset::UInt
end
Base.getindex(s::ShiftedPos, args...) = satadd(s.data[args...], oftype(s.data[args...], s.offset))
Base.length(s::ShiftedPos) = length(s.data)

immutable ShiftedNeg{Q}
	data::Q
	offset::UInt
end
Base.getindex(s::ShiftedNeg, args...) = satsub(s.data[args...], oftype(s.data[args...], s.offset))
Base.length(s::ShiftedNeg) = length(s.data)

immutable Result{P, Q}
	point::P
	r_sq::UInt
	bbox_hi::ShiftedPos{Q}
	bbox_lo::ShiftedNeg{Q}
	Result(point::P) = new(point, typemax(UInt))
	function Result(point::P, r_sq, q::Q)
		r = ceil(Integer, sqrt(r_sq))
		new(point, r_sq, ShiftedPos{Q}(q, r), ShiftedNeg{Q}(q, r))
	end
end

# Euclidean distance, though any p-norm will do.
function sqdist(p, q)
	@assert length(p) == length(q)
	@assert length(q) > 0

	d_sq::UInt = 0
	for i in 1:length(p)
		# d_sq += UInt((p[i] - q[i])^2)
		diff = UInt(p[i] < q[i] ? q[i] - p[i] : p[i] - q[i])
		d_sq = satadd(d_sq, satmul(diff, diff))
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

	# The power-of-two and size of the quadtree-aligned
	# bounding box that most snugly encloses p1 and p2
	power = xor == 0 ? 1 : 1 + exponent(float(xor))
	size = one(Uint) << power

	# Calculate the squared distance from q to the box.
	# The return value is a float for efficiency;
	# it will be multiplied by another float upon return.
	d_sq = 0.0
	for i in 1:length(q)
		# Compute the coordinates of the bounding box
		bbox_lo = (p1[i] >> power) << power
		bbox_hi = bbox_lo + size

		# Accumulate squared distance from the box
		if q[i] < bbox_lo
			d_sq += (Uint(bbox_lo) - q[i])^2
		elseif q[i] > bbox_hi
			d_sq += (q[i] - bbox_hi)^2
		end
	end
	d_sq
end

function nearest{P, Q, T <: Unsigned}(arr::Array{P}, q::Q, lo::T, hi::T, R::Result{P, Q}, ε::Float64)
	# Return early if the range is empty
	lo > hi && return R

	# Calculate the midpoint of the range, avoiding midpoint overflow
	mid = (lo + hi) >>> 1

	# Compute the distance from the probe point to the query point,
	# and update the result if it's closer than our best match so far
	r_sq = sqdist(arr[mid], q)
	r_sq < R.r_sq && (R = Result{P, Q}(arr[mid], r_sq, q))

	# Return early if the range is only one element wide or if the
	# bounding box containing the range is outside of the search radius
	if lo == hi || sqdist_to_quadtree_box(q, arr[lo], arr[hi]) * (1.0 + ε)^2 >= R.r_sq
		return R
	end

	# Recurse. Unlike binary search, we occasionally recurse into
	# both halves of the array when we can't guarantee that the nearest
	# point lies outside the other half
	if shuffless(q, arr[mid])
		R = nearest(arr, q, lo, mid - 1, R, ε)
		shuffmore(R.bbox_hi, arr[mid]) && (R = nearest(arr, q, mid + 1, hi, R, ε))
	else
		R = nearest(arr, q, mid + 1, hi, R, ε)
		shuffless(R.bbox_lo, arr[mid]) && (R = nearest(arr, q, lo, mid - 1, R, ε))
	end

	R
end

function nearest{P, Q}(arr::Array{P}, q::Q, ε=0.0)
	@assert length(arr) > 0 "Searching for the nearest in an empty array"
	nearest(arr, q, UInt(1), UInt(length(arr)), Result{P, Q}(arr[1]), ε).point
end

end # module
