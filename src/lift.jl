# ================================
#   CAKE WAVELETS LIFT
# ================================

# Cut-off function in the frequency plane
mn(ρ::Real; t=1, N = 12) = exp(-ρ^2/t)*sum([(ρ^2/t)^k/factorial(k) for k=0:N])

# B-spline for the radial part
function b_spline(x::T,t1::T,t2::T,t3::T,t4::T) where T<:Real
    if x >= t1 && x<t2
        return (x-t1)^2/(t3-t1)/(t2-t1)
    elseif x >= t2 && x<t3
        return (x-t1)*(t3-x)/(t3-t1)/(t3-t2)+(t4-x)*(x-t2)/(t4-t2)/(t3-t2)
    elseif x >= t3 && x<t4
        return (t4-x)^2/(t4-t2)/(t4-t3)
    else
        return 0
    end
end
b_spline(x) = b_spline(x,-1.,-.5,.5,1.)

# Fourier part of the cake
function cake(ϕ::Real, θs::Int; m::Int = 100, bw::Real = 5, oriented::Bool = false, args...)
    # ϕ = orientation of the cake
    # m = the cake wavelet will be an m x m matrix
    # θs = angular resolution
    # bw = bandwith at which to cut
    @assert m%2 == 0
    
    x(i) = 2bw*(i-1)/m - bw
    ρ(i,j) = sqrt(x(i)^2 +x(j)^2)
    θ(i,j) = atan(x(j),x(i))
    cake(x) = [mn(ρ(i,j))*b_spline((mod(θ(i,j)-x ,2π)- π/2)*θs/(2*pi)) for i in 1:m, j in 1:m]

    if oriented
        cake(ϕ)
    else
        cake(ϕ) + cake(ϕ+π)
    end
end

function lift(img::Array{T,2}, θs::Int = 30 ; m::Int = 100, args...) where T<:Real
    # θs = angular resolution of the cake wvlts
    # m = resolution of the kernel, if 0 it is chosen as the size of the image 
    # (to have perfect reconstruction) 
    @assert(m%2 == 0)
    @assert(size(img, 1)%2 == size(img, 2) %2 == 0)

    m == 0 ? m = minimum(size(img)) : nothing

    lift = zeros(Float64, (size(img)...,θs))
    fimg = img |> fft |> fftshift
    bnd1 = Int(ceil((size(img,1)-m)/2))
    bnd2 = Int(ceil((size(img,2)-m)/2))
    NORM = 1.3333333333333384 # NORM = maximum(cake)
    for k in 1:θs
        flift = zeros(Complex{Float64},size(img))
        cake_wavelet = cake((k-1)*π/θs-π/2, θs; m = m, args...)/NORM
        for r in (bnd1+1):(bnd1 + m), h in (bnd2+1):(bnd2 + m)
            @inbounds  flift[r,h] = cake_wavelet[r-bnd1,h-bnd2]*fimg[r,h]
        end
        @inbounds lift[:,:,k] = flift |> ifftshift |> ifft |> real
    end
    return lift
end
