
module HomogeneousCoordinates

	using Interpolations
	using StaticArrays

	export transform

	"""
		Not sure if this is true for all homogeneous coordinates, but this requires
		the matrix [4,4] element to be equal to one.
	"""
	function transform(a::AbstractArray{<: Real, 3}, b_size::NTuple{N, Int64}, Ma::AbstractMatrix{<: Real}, Mb::AbstractMatrix{<: Real}) where N
		# TODO generated for any dimension
		
		(rows, columns, slices) = size(a)
		

		b = zeros(b_size)
		a_itp = extrapolate(
			interpolate(
				a,
				BSpline(Linear())
			),
			0
		)

		Mba = (inv(Ma) * Mb)[1:3, :]
		for i in CartesianIndices(b_size)
			b_vec = SVector{4, Float64}(i[1], i[2], i[3], 1)
			a_vec = Mba * b_vec
			#a_vec = a_vec ./ a_vec[4] 
			b[i] = a_itp(a_vec...)
		end
		return b
	end

	# This is not yet useful, it can be useful if you only want to change the coordinate system,
	# but not into the system of another scan. In that case the target array size is not defined,
	# and it needs to be determined by transforming the corners. Like that, no empty space is stored.
	# Such a function must then also return the actual indices of the zero index point of the resulting array.
	#b_size = let
	#	corners = @SMatrix [
	#		0 rows	0		rows	0		rows	0		rows;
	#		0 0		columns columns	0		0		columns	columns;
	#		0 0		0		0		slices	slices	slices	slices;
	#		1 1		1		1		1		1		1		1
	#	]
	#	transformed_corners = (inv(Mb) * Ma)[1:3, :] * corners
	#	b_size = map(
	#		limits -> ceil(Int, limits[2] - limits[1]),
	#		dropdims(extrema(transformed_corners; dims=2); dims=2)
	#	)
	#end

end

