using MPIReco, CairoMakie, HDF5
using Statistics
using SparseArrays

include("chebReco/chebReco.jl")

bSF = MPIFile(joinpath(mdfdatadir, "calibrations/13.mdf"), handleSubPeriodsAsFrames=true)
N_ = calibSize(bSF)
N = (N_[1],N_[2],1)  
slice = 2
tfCorrection = false
bgCorrection = true
SMeas = reshape(MPIFiles.getSystemMatrix(bSF; tfCorrection, bgCorrection), N_..., 817, 3)[:,:,slice,:,1:3]

freqs = filterFrequencies(bSF, recChannels = 1:2, minFreq=60e3, maxFreq=450e3)

b_bg = MPIFile(joinpath(mdfdatadir, "measurements/20230613_150948_model based fluid/1.mdf"),handleSubPeriodsAsFrames=true) 
u_bg = getMeasurementsFD(b_bg,numAverages=20, tfCorrection=false)


bs = [ MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/2.mdf"),handleSubPeriodsAsFrames=true),
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/3.mdf"),handleSubPeriodsAsFrames=true), 
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/4.mdf"),handleSubPeriodsAsFrames=true),
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/5.mdf"),handleSubPeriodsAsFrames=true), 
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/6.mdf"),handleSubPeriodsAsFrames=true),
       MPIFile(joinpath(mdfdatadir,"measurements/20230613_150948_model based fluid/7.mdf"),handleSubPeriodsAsFrames=true),
]

spatialScaling = 1.2 #1.2
N = (17,15,1)

p = Dict{Symbol,Any}()
p[:samplingRate] = 2.5e6
p[:amplitude] = Tuple(dfStrength(bSF))
p[:dividers] = (dfDivider(bSF)[1],dfDivider(bSF)[2],1)
p[:gradient] = Tuple(vec(acqGradientDiag(bSF))) .* (-spatialScaling) 
p[:DCore] = 20.0e-9 #18.4e-9 #18.402152475061616*1e-9
p[:kAnis] = 0*2000 #*1800 #900#1800 #0*1835.234880285756 #1835.234880285756
p[:kAnisγ] = 2.0
p[:grid] = RegularGridPositions(collect(N), [0.034, 0.034, 0.006], calibFovCenter(bSF))
#calibFov(bSF), calibFovCenter(bSF))
p[:ensembleAlg] = EnsembleThreads()
p[:model] = EquilibriumModel()
smEq = calcSM(p)
p[:model] = EquilibriumAnisModel()
#smEqAnis = calcSM(p)

KernelUU_aniso, KernelUT_aniso = DCR2_generate_aniso_kernels(p, (21,21))
KernUU_aniso_x = -KernelUU_aniso[:,:,:,1] * 5
KernUU_aniso_y = -KernelUU_aniso[:,:,:,2] * 5
KernUT_aniso_x = -KernelUT_aniso[1:end-1,:,:,1] * 5
KernUT_aniso_y = -permutedims(KernelUT_aniso[1:end-1,:,:,1],(2,1,3)) * 5

#ffn = "../EquilibriumModelPaper/directChebyshevReconstruction/creation_kernel_sys_mat/kernel_anis.h5"
#KernUU_aniso_x_ = -h5read(ffn, "/kernelx") * 4e12
#KernUU_aniso_y_ = -h5read(ffn, "/kernely") * 4e12

tf = -estimateTransferFunction(SMeas, smEq)
#tf = -estimateTransferFunction(SMeas, smEq, tfMeas, shift, phi)
#tf = -estimateTransferFunction(SMeas, smEqAnis)


fig = Figure( size = (1000, 1000*6/length(bs)), figure_padding = 1 )#, fontsize = 16)
#crange = [nothing, nothing,(0,2e-10)]
crange = [nothing, nothing,(0,1.5e7),nothing, nothing,(0,1.5e7)]


for (l,b) in enumerate(bs)

  u = getMeasurementsFD(b,numAverages=100, tfCorrection=false) .- u_bg

  # simulate voltages
  f1 = u[:,1,1,1]
  f2 = u[:,2,1,1]

  fD = 1

  Nb = Int.(p[:dividers][1:2] ./ gcd(p[:dividers][1:2]...))[1] #17



  n_reko = 21;

  f1_ = f1 ./ (tf[:,1].+eps()) #./ (reshape(2*pi*im.*(0:(817-1)),:).+eps())
  f2_ = f2 ./ (tf[:,2].+eps()) #./ (reshape(2*pi*im.*(0:(817-1)),:).+eps())

  ### UU

  c = DCR2(cat(f1_,f2_,dims=2), freqs, [n_reko,n_reko], Nb, fD, "cos", "UU")
  c1 = c[:,:,1]; c2 = c[:,:,2];

  global KernUU_x, KernUU_y = DCR2_generate_kernels(p, (n_reko,n_reko), "UU" )

  # perform deconvolution
  reco_langevin = DCR2_deconvolve_l2(real(c1), real(c2), KernUU_x, KernUU_y, 10^-1.5, p[:gradient][1], p[:gradient][2]);

  if l == 1
    global crange[3] = extrema(real(reco_langevin))
  end

  for (q,reco) in enumerate([c1, c2, reco_langevin])
    ax = CairoMakie.Axis(fig[q, l], ylabel="", xticklabelsvisible=false, xticksvisible=false)
    crange_ = crange[q] == nothing ? extrema(real.(reco)) : crange[q]
    heatmap!(ax, reverse(real.(reco),dims=2), colormap=:inferno,colorrange=crange_)
    hidedecorations!(ax, grid=false, label=false)
    tightlimits!(ax)
  end 

  reko_test_aniso_dir = DCR2_deconvolve_aniso_l2(real(c1), real(c2), 
     KernUU_aniso_x, KernUU_aniso_y, 10^-1.0, 1,1, "direct");
  #reko_test_aniso_dir = DCR2_deconvolve_aniso_l2(real(c2), real(c2), 
  #   KernUU_aniso_y, KernUU_aniso_y, 10^-1.0, 1,1, "direct");
  
  ax = CairoMakie.Axis(fig[4, l], ylabel="", xticklabelsvisible=false, xticksvisible=false)
  heatmap!(ax, reverse(real.(reko_test_aniso_dir),dims=2), colormap=:inferno)#,colorrange=crange_)
  hidedecorations!(ax, grid=false, label=false)
  tightlimits!(ax)  

  ### UT

  #c1ut, c2ut = DCR2(f1_, f2_, freq_comp, n_reko , n_reko , Nb, fD, "cos", "UT");
  cut = DCR2(cat(f1_,f2_,dims=2), freqs, [n_reko,n_reko], Nb, fD, "cos", "UT")
  c1ut = cut[:,:,1]; c2ut = cut[:,:,2];

  global KernUT_x, KernUT_y = DCR2_generate_kernels(p, (n_reko,n_reko), "UT" ) 

  # perform deconvolution
  reco_ut_langevin = DCR2_deconvolve_l2(real(c1ut), real(c2ut), KernUT_x, KernUT_y, 10^-0.7, p[:gradient][1], p[:gradient][2])*2;

  if l == 1
    global crange[6] = extrema(real(reco_ut_langevin))
  end

  for (q,reco) in enumerate([c1ut, c2ut, reco_ut_langevin])
    ax = CairoMakie.Axis(fig[q+4, l], ylabel="", xticklabelsvisible=false, xticksvisible=false)

    crange_ = crange[q+3] == nothing ? extrema(real.(reco)) : crange[q+3]
    heatmap!(ax, reverse(real.(reco),dims=2), colormap=:inferno,colorrange=crange_)
    hidedecorations!(ax, grid=false, label=false)
    tightlimits!(ax)
  end  

  reko_test_aniso_dir = DCR2_deconvolve_aniso_l2(real(c1ut), real(c2ut), 
     KernUT_aniso_x, KernUT_aniso_y, 10^0.5, 1,1, "direct");
  
  ax = CairoMakie.Axis(fig[8, l], ylabel="", xticklabelsvisible=false, xticksvisible=false)
  heatmap!(ax, reverse(real.(reko_test_aniso_dir),dims=2), colormap=:inferno)#,colorrange=crange_)
  hidedecorations!(ax, grid=false, label=false)
  tightlimits!(ax)  

end


strDataset = ["Snail", "Resolution 1", "Resolution 2", "Resolution 3", "Ice Cream", "Dot"]
for l=1:length(strDataset)
  Label(fig[1, l, Top()], strDataset[l], valign = :bottom, font = :bold, #fontsize = 24,
     padding = (0, 0, 5, 0))
end

strModels = ["c1 UU", "c2 UU", "Deconv UU Lang", "Deconv UU EqAniso", "c1 UT", "c2 UT", "Deconv UT Lang", "Deconv UT EqAniso"]

for l=1:length(strModels)
  Label(fig[l, 1, Left()], strModels[l], halign = :left, font = :bold, #fontsize = 24,
    padding = (0, 0, 0, 0), rotation=π/2)
end

rowgap!(fig.layout, 10) 
colgap!(fig.layout, 10)

save("img/recoCheb.png", fig)


fig_,ax,pl = heatmap( reverse(KernUU_x,dims=2), colormap=:inferno, axis=(;title="Kern1"))

heatmap(fig_[1,2], reverse(KernUU_y,dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig_[1,3], reverse(KernUT_x,dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig_[1,4], reverse(KernUT_y,dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig_[2,1], reverse(KernUU_aniso_x[:,:,20],dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig_[2,2], reverse(KernUU_aniso_y[:,:,20],dims=2), colormap=:inferno, axis=(;title="Kern2"))

save("img/kernel.png", fig_)

fig_,ax,pl = heatmap( reverse(KernUU_aniso_x[:,:,1],dims=2), colormap=:inferno, axis=(;title="Kern1"))

heatmap(fig_[1,2], reverse(KernUU_aniso_y[:,:,1],dims=2), colormap=:inferno, axis=(;title="Kern2"))

heatmap(fig_[1,3], reverse(KernUU_aniso_x[:,:,20],dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig_[1,4], reverse(KernUU_aniso_y[:,:,20],dims=2), colormap=:inferno, axis=(;title="Kern2"))

heatmap(fig_[2,1], reverse(KernUU_aniso_x[:,:,40],dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig_[2,2], reverse(KernUU_aniso_y[:,:,40],dims=2), colormap=:inferno, axis=(;title="Kern2"))

heatmap(fig_[2,3], reverse(KernUU_aniso_x[:,:,60],dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig_[2,4], reverse(KernUU_aniso_y[:,:,60],dims=2), colormap=:inferno, axis=(;title="Kern2"))


#heatmap(fig[1,3], reverse(KernUT_x,dims=2), colormap=:inferno, axis=(;title="Kern2"))
#heatmap(fig[1,4], reverse(KernUT_y,dims=2), colormap=:inferno, axis=(;title="Kern2"))

save("img/kernel_aniso.png", fig_)


#=fig,ax,pl = heatmap( reverse(kernels[:,:,1],dims=2), colormap=:inferno, axis=(;title="Kern1"))

heatmap(fig[1,2], reverse(kernels[:,:,10],dims=2), colormap=:inferno, axis=(;title="Kern2"))

heatmap(fig[1,3], reverse(kernels[:,:,20],dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig[1,4], reverse(kernels[:,:,30],dims=2), colormap=:inferno, axis=(;title="Kern2"))

heatmap(fig[2,1], reverse(kernels[:,:,60],dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig[2,2], reverse(kernels[:,:,90],dims=2), colormap=:inferno, axis=(;title="Kern2"))

heatmap(fig[2,3], reverse(kernels[:,:,120],dims=2), colormap=:inferno, axis=(;title="Kern2"))
heatmap(fig[2,4], reverse(kernels[:,:,140],dims=2), colormap=:inferno, axis=(;title="Kern2"))


#heatmap(fig[1,3], reverse(KernUT_x,dims=2), colormap=:inferno, axis=(;title="Kern2"))
#heatmap(fig[1,4], reverse(KernUT_y,dims=2), colormap=:inferno, axis=(;title="Kern2"))

save("img/kernel_aniso2.png", fig)=#

fig