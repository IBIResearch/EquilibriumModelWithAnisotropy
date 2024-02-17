using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MNPDynamics
using CairoMakie, Printf


fb = 2.5e6 # Sampling frequency
fx = fb / 102  # Excitation frequency

samplingMultiplier = 2                  # sampling rate = samplingMultiplier*fb
tLength = samplingMultiplier*102        # length of time vector
tMax = (102-1/samplingMultiplier) / fb  # maximum evaluation time in seconds
t = range(0,stop=tMax,length=tLength)   # time intervall

# Parameters
p = Dict{Symbol,Any}()

p[:α] = 0.1               # damping coefficient
p[:N] = 20                # maximum spherical harmonics index to be considered
p[:relaxation] = NEEL     # relaxation mode
p[:reltol] = 1e-6         # relative tolerance
p[:abstol] = 1e-8         # absolute tolerance
p[:tWarmup] = 2 / fb      # warmup time
p[:derivative] = true     # calculate directly the derivative
p[:order] = 120           # high order for the reference simulation
p[:model] = EquilibriumAnisModel()
dir = [1;0;0]             # colinear easy axis

# Magnetic field for simulation 
const amplitude = 0.012
B =  t -> amplitude*[cospi(2*fx*t), 0, 0]

Ds = range(15e-9, 25e-9, length=29)  # particle core diameter
kAnis = range(0,10000, length=29)    # anisotropy constant
orders = 1:50                        # truncation indices

function calcErr(A,B)
  D = mean(abs.(A.-B), dims=1) ./ maximum(abs.(B), dims=1)
  return D[1]
end

if true # can be set to false to tweak the plot when running the script multiple times
  errors = zeros(length(Ds), length(kAnis), length(orders))
  for d=1:length(Ds)
    for k=1:length(kAnis)
      
      p[:DCore] = Ds[d]
      p[:kAnis] = kAnis[k]*dir
      p[:order] = 200
      mRef_ = simulationMNP(B, t; p...)

      for o = 1:length(orders)
        p[:order] = orders[o]
        mApprox_ = simulationMNP(B, t; p...)

        errors[d,k,o] = calcErr(mApprox_, mRef_)
      end
    end
  end
end

minOrder = [findfirst(d->d<1e-6, errors[d,k,:])  for k=1:length(kAnis), d=1:length(Ds)]


fig = Figure( size = (600, 450), figure_padding = 1 )

ax = CairoMakie.Axis(fig[1, 1], xlabel = rich(rich("D", font=:italic) , " / nm" ), #xlabel=L"D\,/\,\text{nm}",
        ylabel = rich(rich("K", font=:italic), superscript("anis"), " / Jm", superscript("-3") ),  # L"K_\text{anis}\,/\,\text{Jm}^{-3}" ,
        title = rich("Minimal Order for 10", superscript("-6")," Error") #, font=:italic) #L"\text{Minimal Order for } 10^{-6} \text{ Error}",
        ) #titlesize = 20)

hm = CairoMakie.heatmap!(ax, round.(collect(Ds).*1e9, digits=2), kAnis, minOrder,
        colormap = :cividis) 

contour!(ax, round.(collect(Ds).*1e9, digits=2), kAnis, Float64.(minOrder);
       color=:white, linewidth=1.85, levels=[15, 20, 30],
       labels = true, labelsize = 14, labelformatter=(x) -> @sprintf("%d", x))

fig[1,2] = Colorbar(fig, hm,  label = "Order", minorticksvisible = true, minorticks = IntervalsBetween(9))

resize_to_layout!(fig)

save("img/order.png", fig, px_per_unit = 2, pt_per_unit = 0.25)
save("img/order.pdf", fig, px_per_unit = 2, pt_per_unit = 0.25)

fig