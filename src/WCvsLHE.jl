module WCvsLHE

using Images, ImageFiltering
using FFTW
import OffsetArrays

const Lift{T} = Array{T,3}
const Kern{T,n} = OffsetArrays.OffsetArray{T,n,Array{T,n}}

include("lift.jl")

export wc, lhe

normalize(x::Array) = typeof(x)(x - ones(size(x))*minimum(x))/(maximum(x)-minimum(x))

function project(F::Lift; normalized::Bool = true, args...) 
    x = dropdims(sum(F, dims = 3), dims = 3)
    normalized ? x = normalize(x) : nothing
    return x
end

# Approximation of sigmoid for LHE
# Coefficients for the polynomial approximation of the sigmoid
const c = reverse([-7.7456e+00, 
                3.1255e-16, 
                1.5836e+01, 
                -1.8371e-15, 
                -1.1013e+01, 
                4.4531e-16, 
                3.7891e+00, 
                1.2391e-15]);

#Degree of the approximating polynomial
const n = length(c)-1;  

gr(i::Int, img::Array, ω::Kern) = imfilter(img.^i,ω)
a(i::Int, img::Array)=sum( [(-1.)^(j-i+1)*c[j+1]*binomial(j,i)*img.^(j-i) for j in i:n] )
R(img::Array, ω::Kern) = sum( [ a(i,img) .* gr(i,img,ω) for i in 0:n ] )

# Weight functions generation
function adapted_gaussian(σ::Real, I::Array)
    # Generates a Gaussian with variance σ and correct
    # dimensions w.r.t. to I
    if ndims(I) == 2
        var = (σ, σ)
    else
        var = (σ, σ, (size(I, 3)/size(I, 1))*σ)
    end

    # var = tuple([σ for _ in 1:ndims(I)]...)
    Kernel.gaussian(var) |> reflect # Convolution not correlation
end

# Local Mean Average
LMA(σμ::Real, I0::Array{T, 2}) where T = imfilter(I0, adapted_gaussian(σμ, I0))

# Type to save the results
struct Result{n}
    res :: Array{Float64, n}
    iter :: Int
    tol :: Float64
end
Base.getindex(R::Result, etc...)  = getindex(R.res, etc...)
iterations(R::Result) = R.iter
tolerance(R::Result) = R.tol
project(R::Result{3}, etc...) = Result{2}(project(R.res, etc...), iterations(R), tolerance(R))
show(R::Result{2}) = Gray.(R[:,:])

# Actual procedure
function gradient_descent(I0::Array, W, λ::Real, lma::Array; Δt::Real = .15, threshold::Real = .01, max_iter::Int = 50, M::Real = 1, args...)
    # threshold : the iterations stops when two successive terms have a weighted L2 difference smaller than this quantity
    # max_iter: maximum number of iterations
    next(Ik) = (1-(1+λ)Δt)Ik + (lma + W(Ik)/(2*M))Δt + λ*Δt*I0
    proj(x) = sum(x, dims = 3) # it does nothing if x has 2 dims
    diff(a,b) = sqrt( sum(( proj(a-b) ).^2)/sum(proj(a).^2) )

    prec = I0
    cur = next(I0)
    i = 1
    while diff(prec, cur) >= threshold && i < max_iter
        (prec, cur) = (cur, next(cur))
        i += 1
    end

    Result{ndims(I0)}(cur, i, diff(prec, cur))
end

function algo(I0::Array, σμ::Real, σw::Real, λ::Real; algo_type = :planar, model = :wc, θs::Int = 30, α::Int = 5, verbose::Bool = false, args...) 

	if algo_type == :cortical
		F0 = lift(I0, θs; args...)
		lma = lift(LMA(σμ, I0), θs; args...)
	elseif algo_type == :planar
		F0 = copy(I0)
		lma = LMA(σμ, I0)
	else
		error("Type of evolution must be :cortical or :planar")
	end

	if model == :wc
        function σ(x)
            if x < 1/2 - 1/α
                y = 1/2
            elseif x > 1/2 + 1/α
                y = -1/2
            else
                y = (x - 1/2)*α/2
            end
            y
        end

		W(I) = imfilter(σ.(I), adapted_gaussian(σw, F0))
        res = gradient_descent(F0, W, λ, lma; args...)
	elseif model == :lhe
		W2(I) = R(I, adapted_gaussian(σw, F0))
        res = gradient_descent(F0, W2, λ, lma; args...)
    else
		error("Model must be :wc or :lhe")
	end

	if algo_type == :cortical
		res = project(res)
	end	

    verbose ? res : show(res)
end

wc(I0, σμ, σw, λ; args...) = algo(I0, σμ, σw, λ, model = :wc; args...) 
lhe(I0, σμ, σw, λ; args...) = algo(I0, σμ, σw, λ, model = :lhe; args...)

end