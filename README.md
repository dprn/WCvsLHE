# Codes for contrast enhancement via Wilson-Cowan and Local Histogram Equalization

## Usage

In order to import the package, do as usual

```
import Pkg
Pkg.activate(path_to_WCvsLHE_folder)

using WCvsLHE
```

The only exported functions are `wc` (implementing the WC contrast enhancement) and `lhe` (for the Local Histogram Equalization). Both these functions take the following parameters:

- `I0`:	Input image (given as an array of real numbers)
- `σμ`: Standard deviation for the local mean average
- `σw`: Standard deviation for the interaction kernel
- `λ`: Initial data attachment parameter

Moreover, they have the following keyword arguments:

- `M`: Strength of the interaction kernel, which is given as `ν = 1/(2M)` (defaults to `1`)
- `model`: Can be either `:planar` for the non-oriented version or `:cortical` for the oriented one (defaults to `:planar`)
- `θs`: Number of orientations to consider in case of a `:cortical` algorithm (defaults to 30)
- `α`: Sigmoid parameter, used only in `wc` (defaults at `5`)
- `verbose`: If `true` the algorithm returns a `Result` instance, which saves the result, the number of iterations, and the relative tolerance attained. If `false` we only return the final image. (Defaults at `false`)
- `max_iter`: Maximal number of iterations.
- `tolerance`: Tolerance at which to stop.

## Example

```
img = load("myimage.png")
x = Float64.(img)

wc(x, 2, 5, .7, algo_type = :planar)
lhe(x, 2, 5, .7, algo_type = :planar)
wc(x, 2, 5, .7, algo_type = :cortical)
lhe(x, 2, 5, .7, algo_type = :cortical)
```
