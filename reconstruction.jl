using MPIReco, CairoMakie
using Statistics, LinearAlgebra

## SM ##

bSF = MPIFile(joinpath(mdfdatadir, "calibrations/13.mdf"), handleSubPeriodsAsFrames=true)
N_ = calibSize(bSF)
N = (N_[1],N_[2],1)  
slice = 2
tfCorrection = false
bgCorrection = true
SMeas = reshape(MPIFiles.getSystemMatrix(bSF; tfCorrection, bgCorrection), N_..., 817, 3)[:,:,slice,:,1:3]

## Reco params ##

SNRThresh = 0;
λ = 0.1*10^(-0)
iters = 100;
recosize = (17,15)

freqs_ = filterFrequencies(bSF,SNRThresh=SNRThresh,recChannels = 1:2,
                          minFreq=0, maxFreq=1250e3)

freqs = Vector{CartesianIndex{2}}(undef, 0)
similarityThresh = 0.99
for k in freqs_
  B = reshape(ComplexF32.(SMeas)[:,:,:,:], prod(recosize), 817, 3)[:,k]
  A = reshape(ComplexF32.(SCorr[:EqAnis_TFEst])[:,:,:,:], prod(recosize), 817, 3)[:,k]
  err = norm(A-B) / maximum(abs.(B)) / sqrt(length(B))
  if err < similarityThresh
    push!(freqs, k)
  end
end

SMeas_ = transpose(reshape(reshape(ComplexF32.(SMeas)[:,:,:,:], prod(recosize), 817, 3)[:,freqs],prod(recosize),:))

models = [:FP_TFEst, :EqAnis_TFEst, :EqAnisRed_TFEst, :Eq_TFEst]

SCorr_ = [ transpose(reshape(reshape(ComplexF32.(SCorr[m])[:,:,:,:], prod(recosize), 817, 3)[:,freqs],prod(recosize),:)) for m in models]

prepend!(SCorr_, [SMeas_])

## bg

b_bg = MPIFile(joinpath(mdfdatadir, "measurements/20230613_150948_model based fluid/1.mdf"),handleSubPeriodsAsFrames=true) 
u_bg = getMeasurementsFD(b_bg,frequencies=freqs,numAverages=20, tfCorrection=false)

u_bg_ = getMeasurementsFD(b_bg,frequencies=freqs, tfCorrection=false)

bg_std = std(u_bg_, dims=3)
weights = Float32.(mean(abs.(vec(bg_std))) ./ abs.(vec(bg_std)))

## measurements ##

bs = [ MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/2.mdf"),handleSubPeriodsAsFrames=true),
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/3.mdf"),handleSubPeriodsAsFrames=true), 
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/4.mdf"),handleSubPeriodsAsFrames=true),
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/5.mdf"),handleSubPeriodsAsFrames=true), 
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/6.mdf"),handleSubPeriodsAsFrames=true),
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/7.mdf"),handleSubPeriodsAsFrames=true),
]

c_ = Array{Any,2}(undef, length(bs), length(SCorr_))

for (i,S) in enumerate(SCorr_)
  for (l,b) in enumerate(bs)
    u = getMeasurementsFD(b,frequencies=freqs,numAverages=20, frames=1:20,
          spectralLeakageCorrection=false, tfCorrection=false) .- u_bg
    S_ = diagm(weights)*S 
    c = reshape(MPIReco.reconstruction(S_,weights.*u, 
                solver = Kaczmarz,
                reg = [L2Regularization(λ)],
                iterations=iters), recosize[1], recosize[2], 1)
    c_[l,i] = c 
  end
end

fig = Figure( size = (600, 600*(length(SCorr_)+1)/length(bs)), figure_padding = 1 )

fnPhotos = joinpath.(mdfdatadir, ["photos/Snail.jpg", "photos/Resolution1.jpg", 
   "photos/Resolution2.jpg", "photos/Resolution3.jpg", "photos/IceCream.jpg", "photos/Dot.jpg"])

for (l,fn) in enumerate(fnPhotos)
  ax = CairoMakie.Axis(fig[1, l], ylabel="", 
      xticklabelsvisible=false, xticksvisible=false)

  A = imresize( reverse( transpose(load(fn)  ), dims=2), ratio = 0.1)

  CairoMakie.heatmap!(ax, A)  
  hidedecorations!(ax, grid=false, label=false)
  tightlimits!(ax)

  if (l == 1) 
    as_ = 7.6
    fs_ = 14
    arrows!(ax, [size(A,1)*0.05,size(A,1)*0.05], [size(A,2)*0.95,size(A,2)*0.95], 
                [size(A,1)*0.1,0], [0,-size(A,2)*0.1], color=:black,
            arrowsize=as_, lengthscale=3.0, linewidth = 1.5 )

    text!(ax, [size(A,1)*0.02,size(A,1)*(0.4)],[size(A,2)*(0.95-0.45),size(A,2)*(0.87)], text=[rich("x", font=:bold_italic),rich("y", font=:bold_italic)],
            color=:black, fontsize=fs_ )
  end
end


for (i,S) in enumerate(SCorr_)
  for (l,b) in enumerate(bs)
    ax = CairoMakie.Axis(fig[i+1, l], ylabel="", 
        xticklabelsvisible=false, xticksvisible=false)

    A = reverse(squeeze(sum(c_[l,i],dims=3))',dims=2)

    cr = (maximum(A)*0.0,maximum(A))

    CairoMakie.heatmap!(ax, A, colormap=:inferno, colorrange=cr)  
    hidedecorations!(ax, grid=false, label=false)
    tightlimits!(ax)
  end
end

strDataset = ["Snake", "Resolution 1", "Resolution 2", "Resolution 3", "Ice Cream", "Dot"]
for l=1:length(strDataset)
  Label(fig[1, l, Top()], strDataset[l], valign = :bottom, font = :bold, 
       padding = (0, 0, 5, 0))
end

strModels = string.(models)
strModels = ["FP", "EQANIS", "red. EQANIS", "EQ"]

prepend!(strModels, ["Photo", "Measured"])

for l=1:length(strModels)
  Label(fig[l, 1, Left()], strModels[l], halign = :left, font = :bold, 
  padding = (0, 5, 0, 0), rotation=π/2)
end

rowgap!(fig.layout, 5) 
colgap!(fig.layout, 5)

save(joinpath(imgdir,"reco.pdf"), fig)


