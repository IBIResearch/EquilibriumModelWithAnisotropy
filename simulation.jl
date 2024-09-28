using Pkg
Pkg.activate(".")
Pkg.instantiate()

using MNPDynamics
using MPIReco
using CairoMakie
using FFTW
using Serialization
using Statistics
using ImageUtils
using MPIFiles
using LazyArtifacts

const mdfdatadir = joinpath(artifact"MDFStore")
@info "The mdf data is located at $mdfdatadir."

if length(ARGS) == 2
  const datadir = ARGS[1]
  const imgdir = ARGS[2]
else
  const datadir  = mkpath("./data")
  const imgdir  = mkpath("./img")
end

@info "The generated data are located at $datadir."
@info "The generated images are located at $imgdir."

include("utils.jl")
include("systemMatrix.jl")
include("visualization.jl")

datasets = [joinpath(mdfdatadir, "calibrations/13.mdf"), joinpath(mdfdatadir,"calibrations/4.mdf")]
datasetNames = ["SMFluid", "SMImmobilized"]

Ds = [19*1e-9, 19*1e-9]
Ks = [3500.0, 1400.0]
Ksγ = [2.0, 1.0]
dfAmplitudes = [[0.012, 0.012, 0.0], [0.012, 0.012, 0.0]]
scalings = [ [1.15,1.15,1.0], [14.35 / 13.63, 14.71 / 13.63, 1.0]]

α = 1.5
anisotropyAxes = [nothing, (cos(pi/2*α), sin(pi/2*α), 0.0) ]

function simulateAndPlot(dataset, datasetName, D, kAnis, kAnisγ, dfAmplitude, scaling, anisotropyAxis )

  p = Dict{Symbol,Any}()

  tfCorrection = false
  bgCorrection = true
  bSF = MPIFile(dataset, handleSubPeriodsAsFrames=true)
  N_ = calibSize(bSF)
  N = (N_[1],N_[2],1)  
  slice = N_[3] == 1 ? 1 : 2
  p[:anisotropyAxis] = anisotropyAxis

  SMeas = reshape(MPIFiles.getSystemMatrix(bSF; tfCorrection, bgCorrection), N_..., 817, 3)[:,:,slice,:,1:3]


  # Parameters
  p[:DCore] = D  # particle diameter in nm
  p[:α] = 0.1                # damping coefficient
  p[:kAnis] = kAnis         # anisotropy constant
  p[:kAnisγ] = kAnisγ
  p[:N] = 20                 # maximum spherical harmonics index to be considered
  p[:relaxation] = NEEL      # relaxation mode
  p[:model] = FokkerPlanckModel()

  p[:reltol] = 1e-4          # relative tolerance
  p[:abstol] = 1e-6          # absolute tolerance
  p[:tWarmup] = 0.00002      # warmup time
  p[:derivative] = false
  p[:solver] = :FBDF         # Use more stable solver
  p[:samplingRate] = 2.5e6

  p[:amplitude] = dfAmplitude
  p[:dividers] = (dfDivider(bSF)[1],dfDivider(bSF)[2],1)
  p[:grid] = RegularGridPositions(collect(N), calibFov(bSF), calibFovCenter(bSF))
  p[:gradient] = Tuple(vec(acqGradientDiag(bSF))) .* (-scaling) 
  p[:ensembleAlg] = EnsembleThreads()

  tfMeas = rxTransferFunction(bSF)

  filenameSM = joinpath(datadir, "$(datasetName).bin")

  if isfile(filenameSM) && true #false
    sms = deserialize(filenameSM)
  else
    p[:model] = EquilibriumModel()
    smEq = calcSM(p)
    p[:model] = EquilibriumAnisModel()
    smEqAnis = calcSM(p)
    p[:model] = FokkerPlanckModel()
    smFP = calcSM(p)

    smsEqAnisRed = calcSM(p, chebyshev=true)

    sms = Dict{Symbol,Any}()
    sms[:Eq] = smEq
    sms[:EqAnis] = smEqAnis
    sms[:EqAnisRed] = smsEqAnisRed
    sms[:FP] = smFP
    serialize(filenameSM, sms)
  end

  phi = 0.0
  shift = 4.0
  SCorr = Dict{Symbol,Any}(
    :Eq_TFEst => correctTransferFunction(SMeas, sms[:Eq]), # estim
    :Eq_TFMeas => correctTransferFunction(SMeas, sms[:Eq], tfMeas, shift, phi), # measured
    :EqAnis_TFEst => correctTransferFunction(SMeas, sms[:EqAnis]), # estim
    :EqAnis_TFMeas => correctTransferFunction(SMeas, sms[:EqAnis], tfMeas, shift, phi), # measured
    :EqAnisRed_TFEst => correctTransferFunction(SMeas, sms[:EqAnisRed]), # estim
    :EqAnisRed_TFMeas => correctTransferFunction(SMeas, sms[:EqAnisRed], tfMeas, shift, phi), # measured
    :FP_TFEst => correctTransferFunction(SMeas, sms[:FP]), # estim
    :FP_TFMeas => correctTransferFunction(SMeas, sms[:FP], tfMeas, shift, phi) # measured
  )

  plotSMs(SMeas, SCorr, bSF, 1, datasetName, anisotropyAxis; imgDir=imgdir)
  plotSMs(SMeas, SCorr, bSF, 2, datasetName, anisotropyAxis; imgDir=imgdir)
  return SCorr
end

SCorr = simulateAndPlot(datasets[1], datasetNames[1], Ds[1], Ks[1], Ksγ[1], 
                dfAmplitudes[1], scalings[1], anisotropyAxes[1] )

SCorrImmobilized = simulateAndPlot(datasets[2], datasetNames[2], Ds[2], Ks[2], Ksγ[2], 
                dfAmplitudes[2], scalings[2], anisotropyAxes[2] )
