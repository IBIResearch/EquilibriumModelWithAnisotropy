using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MNPDynamics, MPIFiles
using CairoMakie
using Serialization
using Statistics
using BenchmarkTools, Printf

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5

include("utils.jl")

p = Dict{Symbol,Any}()

tfCorrection = false
bgCorrection = true

bSF = MPIFile("data/oriented/SFdeg45.mdf") #,handleSubPeriodsAsFrames=true)
slice = 1
α = 1.5
p[:anisotropyAxis] = (cos(pi/2*α), sin(pi/2*α), 0.0)

N_ = calibSize(bSF)
N = (N_[1],N_[2],1)
@info N

# Parameters
p[:α] = 0.1                # damping coefficient
p[:N] = 20                 # maximum spherical harmonics index to be considered
p[:relaxation] = NEEL      # relaxation mode
p[:model] = FokkerPlanckModel()

p[:kAnisγ] = 1.0 
spatialScaling = 1.0 


p[:reltol] = 2e-4          # relative tolerance
p[:abstol] = 1e-5          # absolute tolerance
p[:tWarmup] = 0.00002      # warmup time
p[:derivative] = false
p[:solver] = :FBDF        
p[:samplingRate] = 2.5e6

p[:amplitude] = Tuple(dfStrength(bSF))
p[:dividers] = (dfDivider(bSF)[1],dfDivider(bSF)[2],1)
p[:grid] = RegularGridPositions(collect(N), calibFov(bSF), calibFovCenter(bSF))
p[:gradient] = Tuple(vec(acqGradientDiag(bSF))) .* (-spatialScaling) 

p[:ensembleAlg] = EnsembleThreads()
#p[:ensembleAlg] = EnsembleSerial()
#p[:ensembleAlg] = EnsembleDistributed()

p[:order] = 45

filenameSM = "data/SMAccuracy.bin"
filenameTimes = "data/SMTimes.bin"

Ds = range(15e-9, 25e-9, length=9)
kAnis = range(0,10000, length=9)

if isfile(filenameSM) && true
  sms = deserialize(filenameSM)
  times = deserialize(filenameTimes)
else
  sms = zeros(ComplexF32, N_[1], N_[2],  817, 3, length(Ds), length(kAnis), 3)
  times = zeros(length(Ds), length(kAnis), 3)

  for d=1:length(Ds)
    for k=1:length(kAnis)
      @info "D=$(Ds[d])  k=$(kAnis[k]) "

      p[:DCore] = Ds[d]
      p[:kAnis] =  kAnis[k]

      p[:model] = FokkerPlanckModel()
      times[d,k,1] = @elapsed sms[:,:,:,:,d,k,1] = calcSM(p)
      p[:model] = EquilibriumAnisModel()
      sms[:,:,:,:,d,k,2] = calcSM(p)
      times[d,k,2] = @belapsed calcSM(p)
      times[d,k,3] = @elapsed sms[:,:,:,:,d,k,3] = calcSM(p, chebyshev=true)

      @info times
      GC.gc()
    end
  end

  serialize(filenameSM, sms)
  serialize(filenameTimes, times)
end

smscorr = copy(sms)

for d=1:length(Ds)
  for k=1:length(kAnis)
    for l=2:3
      a = estimateTransferFunction(sms[:,:,:,1:2,d,k,1], sms[:,:,:,1:2,d,k,l])
      smscorr[:,:,:,1:2,d,k,l] = reshape(a,1,1,size(a,1),size(a,2))  .* sms[:,:,:,1:2,d,k,l]
    end
  end
end


smstd = irfft(sms, 1632, 3);
smstdcorr = irfft(smscorr, 1632, 3);


errors = zeros(length(Ds), length(kAnis), 2)
errorscorr = zeros(length(Ds), length(kAnis), 2)


function calcErr(A,B)
  D = mean(abs.(A.-B), dims=3) ./ maximum(abs.(B), dims=3)
  return maximum(D)
end

for d=1:length(Ds)
  for k=1:length(kAnis)
    errors[d,k,1] = calcErr(smstd[:,:,:,1:2,d,k,2], smstd[:,:,:,1:2,d,k,1])
    errors[d,k,2] = calcErr(smstd[:,:,:,1:2,d,k,3], smstd[:,:,:,1:2,d,k,1])
    errorscorr[d,k,1] = calcErr(smstdcorr[:,:,:,1:2,d,k,2], smstd[:,:,:,1:2,d,k,1])
    errorscorr[d,k,2] = calcErr(smstdcorr[:,:,:,1:2,d,k,3], smstd[:,:,:,1:2,d,k,1])
  end
end


fig = Figure( resolution = (1300, 232), figure_padding = 1 )#, fontsize = 16)

str = ["Error: EQANIS", "Error: red. EQANIS"]

crange = (1e-4,1.0)

for l=1:2

  ax = CairoMakie.Axis(fig[1, l], xlabel = rich(rich("D", font=:italic) , " / nm" ),
         ylabel = (l==1) ? rich(rich("K", font=:italic), superscript("anis"), " / Jm", superscript("-3") ) : "",
         title = "$(str[l])", #xticks = round.(collect(Ds).*1e9, digits=2), 
         #yticks = kAnis,
         #yreversed = true,
         yticklabelsvisible=(l==1), yticksvisible=(l==1))

  global hm = CairoMakie.heatmap!(ax, round.(collect(Ds).*1e9, digits=2),
     kAnis,
     errors[:,:,l], colorscale=log10, colorrange = crange )  

  contour!(ax, round.(collect(Ds).*1e9, digits=2), kAnis, errors[:,:,l];
     color=:white, linewidth=1.85, levels=[1e-3, 3e-3, 1e-2, 1e-1],
     labels = true, labelsize = 14, labelformatter=(x) -> @sprintf("%.3f", x))
end


fig[1,3] = Colorbar(fig[1, 2], hm,  label = "Error", minorticksvisible = true, minorticks = IntervalsBetween(9)) #limits = crange,


strTimes = ["Comp. Time: FP", "Comp. Time: EQANIS"]
crangeT = (1e-2,1000)

for l=1:2

  ax = CairoMakie.Axis(fig[1, l+3], xlabel = rich(rich("D", font=:italic) , " / nm" ),
         ylabel=(l==3) ? rich(rich("K", font=:italic), superscript("anis"), " / Jm", superscript("-3") ) : "",
         title = "$(strTimes[l])", 
         yticklabelsvisible=(l==3), yticksvisible=(l==3))

  global hm2 = CairoMakie.heatmap!(ax, colorscale=log10, round.(collect(Ds).*1e9, digits=2),
     kAnis,
     times[:,:,l], colorrange = crangeT ) 
end

fig[1,6] = Colorbar(fig[1, 4], hm2,  label = rich(rich("t", font=:italic), " / s"), minorticksvisible = true, minorticks = IntervalsBetween(9)) #limits = crange,

crangeS = (1e1,1e4)
ax = CairoMakie.Axis(fig[1, 7], xlabel = rich(rich("D", font=:italic) , " / nm" ),
         ylabel=false ? rich(rich("K", font=:italic), superscript("anis"), " / Jm", superscript("-3") ) : "",
         title = "Speedup: EQANIS / FP", 
         yticklabelsvisible=false, yticksvisible=false)

hm3 = CairoMakie.heatmap!(ax, round.(collect(Ds).*1e9, digits=2),
     kAnis, colorscale=log10, colorrange = crangeS,
     times[:,:,1]./times[:,:,2])

fig[1,8] = Colorbar(fig, hm3,  label = "Speedup", minorticksvisible = true, minorticks = IntervalsBetween(9))

save("img/accuracy.pdf", fig, px_per_unit = 2, pt_per_unit = 0.25)