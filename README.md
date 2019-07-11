# Codes for contrast enhancement via Wilson-Cowan and Local Histogram Equalization

## How to reproduce the experiments of the paper

From inside the `experiments` folder, it suffices to run the following command

```sh
julia experiments.jl
```

## Usage

In order to import the package, do as usual

```julia
import Pkg
Pkg.activate(path_to_WCvsLHE_folder)

using WCvsLHE
```

The only exported functions are `wc` (implementing the WC contrast enhancement) and `lhe` (for the Local Histogram Equalization). There is a special type to store parameters, `Params`, which should be used to call both these functions. It is defined as follows:

```julia
struct Params
    σμ :: Float64 # std deviation of local mean average
    σw :: Float64 # std deviation of interaction kernel
    λ :: Float64 # attachment to the data
    M :: Float64 # the interaction weight is ν=1/(2M)
end
```

The functions `wc` and `lhe` take the following parameters:

- `I0`:	Input image
- `p`: Parameters

Moreover, they have the following keyword arguments:

- `model`: Can be either `:planar` for the non-oriented version or `:cortical` for the oriented one (defaults to `:planar`)
- `θs`: Number of orientations to consider in case of a `:cortical` algorithm (defaults to 30)
- `α`: Sigmoid parameter, used only in `wc` (defaults at `5`)
- `verbose`: If `true` the algorithm returns a `Result` instance, which saves the result, the number of iterations, and the relative tolerance attained. If `false` we only return the final image. (Defaults at `false`)
- `max_iter`: Maximal number of iterations.
- `tolerance`: Tolerance at which to stop.

### Example

```julia
using Images
img = load("myimage.png")

p = Params(2, 5, .7, .7)

wc(img, p, algo_type = :planar)
lhe(img, p, algo_type = :planar)
wc(img, p, algo_type = :cortical)
lhe(img, p, algo_type = :cortical)
```

