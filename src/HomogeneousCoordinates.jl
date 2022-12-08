
module HomogeneousCoordinates

	using Base.Cartesian
	using StaticArrays
	using LinearAlgebra
	using Interpolations

	"""
		matrices [N,N] elements must be equal to one.
	
		The shifts in the fourth column of the matrices must be
		computed from the centre of the array, defined as 1 .+ size(a) .÷ 2
		Note:
		This is more useful than a shift from the index origin (1,1,...), because:
			- if the shifts of two arrays are the same, then their origins match.
			- this shift if the shift of a symmetric volume around the origin, unpolluted by the array size

	"""
	# TODO mutating for multi channel
	function transform(
		a::AbstractArray{T, N},
		b_size::NTuple{N, Int64},
		Ma::AbstractMatrix{<: Real},
		Mb::AbstractMatrix{<: Real},
		interp_scheme::Interpolations.InterpolationType,
		extrap_scheme::Union{Int, Interpolations.BoundaryCondition}
	) where {N, T <: Number}
		# TODO generated for any dimension
		
		b = zeros(T, b_size)
		a_itp = extrapolate(
			interpolate(
				a,
				interp_scheme
			),
			extrap_scheme
		)

		# Get joint transformation matrix
		#=
			Compute offset due to one based indexing and origin of the array

			Example, say a and b have the number of elements but b's axis goes in the opposite direction.
			Both shifts are zero.

			a (size = 6):
			1	2	3	4	5	6	(indices)
						[			(origin)
			1	0	0	0	0	0	(values)
					

			and thus,
			b (size = 6)
				6	5	4	3	2	1	(indices)
						]				(origin)
				0	0	0	0	0	0	(values)

			And consequently the one in `a` would not be interpolated into `b`.


			index_a	= Ma^-1 * Mb * ( index_b - [1, ..., 0]  - size(b) ÷ 2 ) +
						+ [1, ..., 0] + size(a) ÷ 2
					= Ma^-1 * Mb * index_b + 
						- Ma^-1 * Mb * ([1, ..., 0] + size(b) ÷ 2 ) +
						+ [1, ..., 0] + size(a) ÷ 2

			This shift can be added to the shift in the transformation matrix.
		=#
		Mba = let Mba
			Mba = (inv(Ma) * Mb)[1:N, :]
			a_shift = collect(Float64, 1. .+ size(a) .÷ 2)
			b_shift = collect(Float64, 1. .+ b_size .÷ 2)
			@views shift = a_shift - Mba[1:N, 1:N] * b_shift
			Mba[1:N, N+1] += shift
			SMatrix{N, N+1}(Mba)
		end

		# Iterate indices of b
		Threads.@threads for i in CartesianIndices(b_size)
			b_vec = SVector{N+1, Float64}(Tuple(i)..., 1)
			a_vec = Mba * b_vec
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


	"""
		not sure if useful
		Columns are basis vectors
		Make all vectors point into a "positive" direction relative to the canonical basis
	"""
	function adjust_lhrh!(e::AbstractMatrix{<: Real})
		i = argmax(abs2.(e); dims=1)
		signs = sign.(e[i])
		e .*= signs
		return
	end

	"""
		R: volume rotation matrix
		δx: voxel size
		Δx: translation of volume
	"""
	function transformation_matrix(R::AbstractMatrix{<: Real}, δx::AbstractVector{<: Real}, Δx::AbstractVector{<: Real})
		n = size(R, 1)
		@assert n == size(R, 2) == length(δx) == length(Δx)
		# Allocate space
		M = Matrix{Float64}(undef, n+1, n+1)
		# Scale
		S = diagm(δx)
		# Rotation
		M[1:n, 1:n] = R * S
		# Translate
		M[1:n, n+1] .= Δx
		# Fill residual elements
		M[n+1, 1:n] .= 0
		M[n+1, n+1] = 1
		return M
	end

	function scale!(M::AbstractMatrix{<: Real}, δx::AbstractVector{<: Real})
		n = size(M, 1)
		@assert n == size(M, 2) == length(δx)+1
		for i = 1:n-1
			M[1:n-1, i] .*= δx[i]
		end
		return
	end

	# TODO: This functionality doesn't fully belong here, also there is a RotationMatrix package in the standard registry
	function R_x(α::Real)
		sine, cosine = sincos(α)
		return @SMatrix [
			1		0	0;
			0	cosine	-sine;
			0	sine	cosine
		]
	end
	function R_y(α::Real)
		sine, cosine = sincos(α)
		return @SMatrix [
			cosine	0	sine;
			0		1	0;
			-sine	0	cosine
		]
	end
	function R_z(α::Real)
		sine, cosine = sincos(α)
		return @SMatrix [
			cosine	-sine	0;
			sine	cosine	0;
			0		0		1;
		]
	end

	"""
		axis in spherical coordinates (ϕ, θ) points along the third axis of a volume,
		α gives rotation of that volume around this axis
		(first rotation around second axis by θ, then rotation around third axis by ϕ)
		Not tested thoroughly
	"""
	@inline volume_rotation(ϕ::Real, θ::Real, α::Real) = R_y(θ) * R_z(ϕ) * R_z(α) * R_z(-ϕ) * R_y(-θ)
	function cosines2spherical(dc::AbstractVector{<: Real})
		@assert length(dc) == 3
		φ = atan(dc[2], dc[1])
		θ = atan(sqrt(dc[1]^2 + dc[2]^2), dc[3])
		return φ, θ
	end

end

