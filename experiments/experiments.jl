import Plots
using Images

import Pkg
Pkg.activate("../../WCvsLHE")
using WCvsLHE

normalize(x) = typeof(x)( (x-minimum(x).*ones(size(x)))*0.7./(maximum(x)-minimum(x)) + 0.15*(ones(size(x))) )

function plot_comparison_middle(A, B,C;args...)  
    mA = floor(size(A,1)/2) |> Int64
    mB = floor(size(B,1)/2) |> Int64
    mC = floor(size(C,1)/2)|> Int64
    Plots.plot([Float64.(A[mA,:]),Float64.(B[mB,:]), Float64.(C[mC,:])]; args...)
end

function test(name, params; lhe3d = false, dest = name)
	println(string("Processing: ", name))
	img = Gray.(load(string("originals/",name,".png")))
	img = Gray.(normalize(Float64.(img)))
	print(string("\tContrast enhancement: WC2D "))
	wc2d = wc(img, params[:WC2D], algo_type = :planar) |> normalize
	print(string("LHE2D "))
	lhe2d = lhe(img, params[:LHE2D], algo_type = :planar) |> normalize
	plot2d = plot_comparison_middle(img,wc2d,lhe2d, ylims=[.15,.85], linewidth=1.5 ,labels = ["original", "wc2D", "lhe2D"], linecolor = ["blue" "red" "green"], linestyle = [:dot :solid :solid])  
	print(string("WC3D "))
	wc3d = wc(img, params[:WC3D], algo_type = :cortical) |> normalize


	if lhe3d
		print(string("LHE3D "))
		lhe3d = lhe(img, params[:LHE3D], algo_type = :cortical) |> normalize
		plot3d = plot_comparison_middle(img,wc3d,lhe3d, ylims=[.15,.85], linewidth=1.5 ,labels = ["original", "wc3D", "lhe3D"], linecolor = ["blue" "red" "green"], linestyle = [:dot :solid :solid])  
		save(string("results/",dest,"_lhe3d.png"), map(clamp01nan,lhe3d))
		save(string("results/",dest,"_p3d.png"), plot3d)
	end

	println("\n\tSaving.")

	save(string("results/",dest,".png"), img)
	save(string("results/",dest,"_wc2d.png"), map(clamp01nan, wc2d))
	save(string("results/",dest,"_lhe2d.png"), map(clamp01nan, lhe2d))
	save(string("results/",dest,"_wc3d.png"), map(clamp01nan, wc3d))
	save(string("results/",dest,"_p2d.png"), plot2d)
end

# Illusion parameters
illusions = Dict(
"white" => Dict(:WC2D => Params(10,20,.7,1.4),
				:LHE2D => Params(10,50,.7,1),
				:WC3D => Params(20,30,.7,1.4),
				:LHE3D => Params(2, 50,.7,1))

, "brightness" => Dict(:WC2D => Params(2,10,.7,1.4),
				:LHE2D => Params(2,10,.7,1),
				:WC3D => Params(2,10,.7,1.4),
				:LHE3D => Params(2, 10,.7,1))

, "checkerboard" => Dict(:WC2D => Params(10,70,.7,1.4),
				:LHE2D => Params(10,70,.7,1),
				:WC3D => Params(10,70,.7,1.4),
				:LHE3D => Params(10,70,.7,1))

, "chevreul" => Dict(:WC2D => Params(2, 5,.7,1),
				:LHE2D => Params(2,10,.7,1),
				:WC3D => Params(2, 40,.5,1),
				:LHE3D => Params(5,7,.7,1))

, "chevreulCanc" => Dict(:WC2D => Params(2, 20,.5,1.4),
				:LHE2D => Params(2, 20,.5,1),
				:WC3D => Params(2, 20,.5,1.4),
				:LHE3D => Params(2, 20,.5,1))

, "dungeon" => Dict(:WC2D => Params(6,10,.7,1.4),
				:LHE2D => Params(5,40,.7,1),
				:WC3D => Params(2,50,.7,1.4),
				:LHE3D => Params(5,50,.7,1))

, "gratings" => Dict(:WC2D => Params(2,6,.7,1),
				:LHE2D => Params(2,6,.7,1),
				:WC3D => Params(2,6,.7,1),
				:LHE3D => Params(2, 6,.7,1))

, "hong_shevell" => Dict(:WC2D => Params(5,20,.7,1),
				:LHE2D => Params(5,.5,.7,1),
				:WC3D => Params(10,30,.7,1),
				:LHE3D => Params(10,30,.7,1))

, "luminance" => Dict(:WC2D => Params(2,6,.7,1),
				:LHE2D => Params(2,6,.7,1),
				:WC3D => Params(2,6,.7,1),
				:LHE3D => Params(2,6,.7,1))

, "poggendorff" => Dict(:WC2D => Params(3,10,.5,1),
				:LHE2D => Params(3,10,.5,1),
				:WC3D => Params(3,10,.5,1),
				:LHE3D => Params(3,10,.5,1))

# , "tilt" => Dict(:WC2D => Params(15,20,.7,1),
# 				:LHE2D => Params(15,20,.7,1),
# 				:WC3D => Params(15,20,.7,1),
# 				:LHE3D => Params(15,20,.7,1))
)

function batch_test(dict; args...)
	for k in keys(dict)
		if k ∉ ["hong_shevell", "tilt", "dungeon"]
			test(k, dict[k]; args...)
		end
	end

	test("hs1", dict["hong_shevell"]; args...)
	test("hs2", dict["hong_shevell"]; args...)

	test("dungeon1", dict["dungeon"]; args...)
	test("dungeon2", dict["dungeon"]; args...)
end

mkpath("results")
batch_test(illusions, lhe3d = false)
