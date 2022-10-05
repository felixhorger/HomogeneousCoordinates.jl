
module HomogeneousCoordinates

	using Interpolations
	using StaticArrays
	using Base.Cartesian

	"""
		Not sure if this is true for all homogeneous coordinates, but this requires
		the matrix [N,N] element to be equal to one.
	"""
	@generated function transform(
		a::AbstractArray{T, 3},
		b_size::NTuple{N, Int64},
		Ma::AbstractMatrix{<: Real},
		Mb::AbstractMatrix{<: Real},
		interp_scheme::Interpolations.InterpolationType,
		extrap_scheme::Union{Int, Interpolations.BoundaryCondition}
	) where {N, T <: Number}
		# TODO generated for any dimension
		
		return quote
			b = zeros(T, b_size)
			a_itp = extrapolate(
				interpolate(
					a,
					interp_scheme
				),
				extrap_scheme
			)

			# Get joint transformation matrix
			Mba = (inv(Ma) * Mb)[1:$N, :]
			# Compute offset due to one based indexing and origin of the array
			let
				b_shift = Vector{Float64}(undef, $N+1)
				b_shift[1:3] = collect(Float64, 1 .+ b_size .÷ 2)
				b_shift[4] = 1
				@views shift = collect(Float64, 1 .+ size(a) .÷ 2) - Mba[1:3, :] * b_shift
				Mba[1:3, $N+1] += shift
			end

			# Iterate indices of b
			for i in CartesianIndices(b_size)
				b_vec = SVector{$N+1, Float64}(Tuple(i)..., 1)
				a_vec = Mba * b_vec
				#a_vec = a_vec ./ a_vec[N]
				b[i] = $(Expr(:call, :a_itp, ntuple(d -> :(a_vec[$d]), N)...))
			end
			return b
		end
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

	function R_y(α::Real)
		sine, cosine = sincos(α)
		return [
			cosine	0	sine;
			0		1	0;
			-sine	0	cosine
		]
	end
	function R_z(α::Real)
		sine, cosine = sincos(α)
		return [
			cosine	-sine	0;
			sine	cosine	0;
			0		0		1;
		]
	end

	"""
		axis in spherical coordinates (ϕ, θ) points along the third axis of a volume,
		α gives rotation of that volume around this axis
		(first rotation around second axis by θ, then rotation around third axis by ϕ)
	"""
	@inline volume_rotation(ϕ::Real, θ::Real, α::Real) = R_z(ϕ) * R_y(θ) * R_z(α)

	function cosines2spherical(dc::AbstractVector{<: Real})
		@assert length(dc) == 3
		φ = atan(dc[2], dc[1])
		θ = atan(sqrt(dc[1]^2 + dc[2]^2), dc[3])
		return φ, θ
	end

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
		# Rotation
		M[1:n, 1:n] = R
		# Scale
		scale!(M, δx)
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

end

