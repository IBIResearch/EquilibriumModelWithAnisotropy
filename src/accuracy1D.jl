using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MNPDynamics
using CairoMakie, Printf

# Excitation frequencies
fb = 2.5e6
fx = fb / 102

samplingMultiplier = 2                  # sampling rate = samplingMultiplier*fb
tLength = samplingMultiplier*102        # length of time vector
tMax = (102-1/samplingMultiplier) / fb  # maximum evaluation time in seconds
t = range(0,stop=tMax,length=tLength);

# Parameters
p = Dict{Symbol,Any}()

p[:α] = 0.1               # damping coefficient
p[:N] = 20                # maximum spherical harmonics index to be considered
p[:relaxation] = NEEL     # relaxation mode
p[:reltol] = 1e-6         # relative tolerance
p[:abstol] = 1e-8         # absolute tolerance
p[:tWarmup] = 2 / fb      # warmup time
p[:derivative] = true
p[:order] = 200

# Magnetic field for simulation 
const amplitude = 0.012
B =  t -> amplitude*[cospi(2*fx*t), 0, 0]

p[:DCore] = 20e-9         # particle diameter in nm
dirs = [[1;0;0],normalize([0.0;1.0;0])]
Ks = range(0, 10000, length=11)
models = [FokkerPlanckModel(),EquilibriumAnisModel()]

fig = Figure( size = (1300, 280), figure_padding = 1 )
strModels = ["FP", "EQANIS"]
strDirs = ["Colinear Easy Axis", "Orthogonal Easy Axis"]

for d=1:length(dirs)
  for m=1:length(models)
    global ax = CairoMakie.Axis(fig[m,d+1], xlabel= (m==2) ? rich(rich("t", font=:italic), " / μs")  : "", 
         ylabel = (d==1) ? rich(rich("∂", subscript("t"), " m", subscript("x"), font=:italic), " / mAm", superscript("-3"), "s", superscript("-1")) : "",
         title = (m==1) ? "$(strDirs[d])" : "",
         xticklabelsvisible=(m==2), xticksvisible=(m==2))

    for k=1:length(Ks)
      p[:kAnis] = Ks[k]*dirs[d]
      p[:model] = models[m]

      m_ = simulationMNP(B, t; p...)

      lines!(ax, t*1e6, m_[:,1]*1e-6, label = "$(round(Int,Ks[k]))", 
        color = RGBAf(0/255,73/255,146/255,1-0.8*k/length(Ks)), linewidth=1.8)

    end
    tightlimits!(ax, Left(), Right())
  end
end


for l=1:length(strModels)
  Label(fig[l, 1], strModels[l], halign = :left, font = :bold, 
  padding = (0, 0, 0, 0), rotation=π/2, tellheight = false)
end

fig[1:2, 4] = Legend(fig, ax, 
  rich(rich("K", font=:italic), superscript("anis"), " / Jm", superscript("-3"), font=:regular ),
  tellheight=true, 
  rowgap=0.5, patchlabelgap=2, margin =(2.0f0, 2.0f0, 2.0f0, 2.0f0),
  padding = 3, titlegap = 3)


Ds = range(15e-9, 25e-9, length=30)
kAnis = range(0, 10000, length=30)
offsets = range(0, 1.0, length=10)


function calcErr(A,B)
  D = mean(abs.(A.-B), dims=1) ./ maximum(abs.(B), dims=1)
  return D[1] 
end

if true # can be set to false to tweak the plot when running the script multiple times
  errors = zeros(length(Ds), length(kAnis), length(offsets))

  for d=1:length(Ds)
    for k=1:length(kAnis)
      p[:DCore] = Ds[d]
      p[:kAnis] = kAnis[k]*dirs[1]

      for o = 1:length(offsets)
        global B = t -> amplitude*[cospi(2*fx*t)+offsets[o], 0, 0]

        p[:model] = models[1]
        mRef_ = simulationMNP(B, t; p...)
        p[:model] = models[2]
        mApprox_ = simulationMNP(B, t; p...)
        errors[d,k,o] = calcErr(mApprox_, mRef_)
      end
    end
  end
end

crange = (1e-4,10.0)

ax = CairoMakie.Axis(fig[1:2, 5], xlabel = rich(rich("D", font=:italic) , " / nm" ),
        ylabel=  rich(rich("K", font=:italic), superscript("anis"), " / Jm", superscript("-3") ),
        title = "Error: Colinear Easy Axis")

hm = CairoMakie.heatmap!(ax, round.(collect(Ds).*1e9, digits=2), kAnis, 
                 squeeze(maximum(errors,dims=3)), 
                 colorscale=log10, colorrange = crange )  

contour!(ax, round.(collect(Ds).*1e9, digits=2), kAnis, squeeze(maximum(errors,dims=3));
       color=:white, linewidth=1.85, levels=[1e-3, 3e-3, 1e-2, 1e-1],
       labels = true, labelsize = 14, labelformatter=(x) -> @sprintf("%.3f", x))

fig[1:2,6] = Colorbar(fig, hm,  label = "Error", minorticksvisible = true, minorticks = IntervalsBetween(9)) #limits = crange,

colsize!(fig.layout, 5, Aspect(1, 2.2))

resize_to_layout!(fig)

save("img/comparison1D.pdf", fig, px_per_unit = 2, pt_per_unit = 0.25)

fig