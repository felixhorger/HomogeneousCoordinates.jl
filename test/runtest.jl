
using Revise
using LinearAlgebra
using Interpolations
import HomogeneousCoordinates

Ma = collect(Float64, I(2))
Mb = copy(Ma)
Mb[1,1] = -1

a = zeros(6)
b = zeros(6)
a[4] = 1

Base.active_repl.options.iocontext[:displaysize] = (100, 80)
c = HomogeneousCoordinates.transform(
	a,
	(6,),
	Ma, Mb,
	Interpolations.BSpline(Interpolations.Linear()), 0
)

@time HomogeneousCoordinates.R_z(0)

