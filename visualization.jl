function plot2DSM(grid, smFT, MX, MY, rec; ref=nothing, showaxis=false, anisotropyAxis=nothing)

  freqidx = Int[]
  for (iy,my) in enumerate(MY)
    for (ix,mx) in enumerate(MX)
      ax = CairoMakie.Axis(grid[ix, iy], ylabel="", 
          xticklabelsvisible=false, xticksvisible=false)

      idx = (mx-1)*16+(my-1)*17+1
      push!(freqidx, idx)

      c = smFT[:, :, idx, rec]

      if ref == nothing
        ref = copy(smFT)
      end

      if (mx == 1 && my == 1) || (mx == 2 && my == 1) || (mx == 1 && my == 2)
        c .= maximum(abs.(ref[:, :, :, rec]))
      end

      amax = maximum(abs.(ref[:, :, idx, rec]))

      A = collect(transpose(ImageUtils.complexColoring(collect((c)); amax))) #

      CairoMakie.heatmap!(ax, A)  
      hidedecorations!(ax, grid=false, label=false)
      tightlimits!(ax)

      if showaxis && (ix == 1) && (iy == 1)
        as_ = 7.6
        fs_ = 14
        Q1_ = size(A,1)
        Q2_ = size(A,2)
        arrows!(ax, [Q1_*0.1,Q1_*0.1], [Q2_*0.95,Q2_*0.95], 
                    [Q1_*0.1,0], [0,-Q2_*0.1], color=:white,
                arrowsize=as_, lengthscale=3.0, linewidth = 1.5 )
    
        text!(ax, [Q1_*0.05,Q1_*(0.45)],[Q2_*(0.95-0.6),Q2_*(0.82)], text=[rich("x", font=:bold_italic),rich("y", font=:bold_italic)],
                color=:white, fontsize=fs_ )
      end

      if anisotropyAxis != nothing && (ix == 1) && (iy == 1)
        Q1_ = size(A,1)
        Q2_ = size(A,2)
        F = 2
        lines!(ax, [1+F,Q1_-F],[Q2_-F,1+F], color=:yellow, linewidth = 1.5)
      end
    end
  end

  rowgap!(grid, 5) 
  colgap!(grid, 5)

  return freqidx
end




function plotSMs(SMeas, SCorr, bSF, rec = 2, filenamePrefix = "SM", anisotropyAxis = nothing)

  MX = 2:2:9
  MY = 2:2:9

  fig = Figure( size = (1300, 1200), figure_padding = 1 )#, fontsize = 16)

  gr1 = GridLayout()
  fig[1,1] = gr1
  gr2 = GridLayout()
  fig[2,1] = gr2
  gr3 = GridLayout()
  fig[2,2] = gr3
  gr4 = GridLayout()
  fig[2,3] = gr4
  gr5 = GridLayout()
  fig[2,4] = gr5  
  gr6 = GridLayout()
  fig[3,1] = gr6
  gr7 = GridLayout()
  fig[3,2] = gr7
  gr8 = GridLayout()
  fig[3,3] = gr8
  gr9 = GridLayout()
  fig[3,4] = gr9

  smEq = SCorr[:Eq_TFEst]
  smEqAnis = SCorr[:EqAnis_TFEst]
  smEqAnisRed = SCorr[:EqAnisRed_TFEst]
  smFP = SCorr[:FP_TFEst]

  # change to measured TFs
  #smEq = SCorr[:Eq_TFMeas]
  #smEqAnis = SCorr[:EqAnis_TFMeas]
  #smEqAnisRed = SCorr[:EqAnisRed_TFMeas]
  #smFP = SCorr[:FP_TFMeas]

  smModels = [smFP, smEqAnis, smEqAnisRed, smEq]

  freqidx = plot2DSM(gr1, SMeas, MX, MY, rec; showaxis=true, anisotropyAxis)
  plot2DSM(gr2, smFP, MX, MY, rec )
  plot2DSM(gr3, smEqAnis, MX, MY, rec )
  plot2DSM(gr4, smEqAnisRed, MX, MY, rec )
  plot2DSM(gr5, smEq, MX, MY, rec )
  plot2DSM(gr6, smFP - SMeas, MX, MY, rec; ref=SMeas*0.5 )
  plot2DSM(gr7, smEqAnis - SMeas, MX, MY, rec; ref=SMeas*0.5 )
  plot2DSM(gr8, smEqAnisRed - SMeas, MX, MY, rec; ref=SMeas*0.5 )
  plot2DSM(gr9, smEq - SMeas, MX, MY, rec; ref=SMeas*0.5 )

  for r in [(1,1),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4)]
    ax_ = CairoMakie.Axis(fig[r[1],r[2]], xlabel=rich("κ", subscript("y"), font=:italic), 
       ylabel = rich("κ", subscript("x"), font=:italic), yreversed = true, 
       ylabelpadding = 2, xlabelpadding = 2, xticks = MX.-1, yticks = MY.-1, 
       limits=(0.1,maximum(MX)-1+0.9,0.1,maximum(MY)-1+0.9))
  end

  Label(fig[1, 1, Top()], "Measured", valign = :bottom, font = :bold, #fontsize = 24,
     padding = (0, 0, 5, 0))
  strDataset = ["FP", "EQANIS", "red. EQANIS", "EQ"]
  for l=1:length(strDataset)
    Label(fig[2, l, Top()], strDataset[l], valign = :bottom, font = :bold, #fontsize = 24,
         padding = (0, 0, 5, 0))
  end

  Label(fig[3, 1, Top()], "Diff: FP - Measured", valign = :bottom, font = :bold, #fontsize = 24,
         padding = (0, 0, 5, 0))
  Label(fig[3, 2, Top()], "Diff: EQANIS - Measured", valign = :bottom, font = :bold, #fontsize = 24,
         padding = (0, 0, 5, 0))
  Label(fig[3, 3, Top()], "Diff: red. EQANIS - Measured", valign = :bottom, font = :bold, #fontsize = 24,
         padding = (0, 0, 5, 0))
  Label(fig[3, 4, Top()], "Diff: EQ - Measured", valign = :bottom, font = :bold, #fontsize = 24,
         padding = (0, 0, 5, 0))

  C = collect(range(0,1,length=100)) *
  transpose(exp.(2*π*im*collect(range(0,1,length=100))))
  axleg, pleg = CairoMakie.heatmap(fig[2,5],
      (complexColoring(C)),
      axis=(ylabel="Phase", xlabel="Amplitude", #title="colorbar", titlefont = :bold,
      #titlesize = 24, 
      yaxisposition = :right,
      xticks=([1,100],["0","1"]),
      yticks=([0.5,50,100.5],["0", "π", "2π"])))
  axleg, pleg = CairoMakie.heatmap(fig[3,5],
      (complexColoring(C)),
      axis=(ylabel="Phase", xlabel="Amplitude", #title="colorbar", titlefont = :bold,
      #titlesize = 24, 
      yaxisposition = :right,
      xticks=([1,100],["0","0.5"]),
      yticks=([0.5,50,100.5],["0", "π", "2π"])))
  colsize!(fig.layout, 5, Aspect(1, 0.1))

  energyMeas = maximum(reshape(abs.(SMeas),size(SMeas,1)*size(SMeas,2),size(SMeas,3),size(SMeas,4)),dims=1)
  energyFP = maximum(reshape(abs.(smFP),size(SMeas,1)*size(SMeas,2),size(SMeas,3),size(SMeas,4)),dims=1)
  energyEq = maximum(reshape(abs.(smEq),size(SMeas,1)*size(SMeas,2),size(SMeas,3),size(SMeas,4)),dims=1)
  energyEqAnis = maximum(reshape(abs.(smEqAnis),size(SMeas,1)*size(SMeas,2),size(SMeas,3),size(SMeas,4)),dims=1)
  energyEqAnisRed = maximum(reshape(abs.(smEqAnisRed),size(SMeas,1)*size(SMeas,2),size(SMeas,3),size(SMeas,4)),dims=1)

  ax = CairoMakie.Axis(fig[1, 2:5], yscale = log10, alignmode = Outside(10,0,-40,0),
                    xlabel = rich(rich("f", font=:italic), " / kHz"),
                    ylabel = "Row Energy")

  freq = rxFrequencies(bSF) ./ 1000
  fidx = 1:(div(length(rxFrequencies(bSF)),3))

  CairoMakie.lines!(ax, freq[fidx], energyMeas[1,fidx,rec], color = RGBf(colors[1]...), label = "Measured", linewidth=2 ) 
  CairoMakie.lines!(ax, freq[fidx], energyFP[1, fidx,rec], color = RGBf(colors[2]...), label = "FP", linewidth=2 ) 
  CairoMakie.lines!(ax, freq[fidx], energyEqAnis[1,fidx,rec], color = RGBf(colors[3]...), label = "EQANIS", linewidth=2 ) 
  CairoMakie.lines!(ax, freq[fidx], energyEqAnisRed[1,fidx,rec], color = RGBf(colors[4]...), label = "red. EQANIS", linewidth=2 ) 
  CairoMakie.lines!(ax, freq[fidx], energyEq[1,fidx,rec], color = RGBf(colors[5]...), label = "EQ", linewidth=2 ) 

  CairoMakie.scatter!(ax, freq[freqidx], energyMeas[1,freqidx,rec], color = RGBf(colors[1]...), 
     markersize = 14, marker = :xcross ) 

  tightlimits!(ax)
  axislegend(ax)

  # Error Plots
  crange = (0.05,0.36)
  maxMixingError = max(9, MX[end])
  for l=1:length(smModels)

    err = relErrorMixingOrder(smModels[l], SMeas, 1:maxMixingError )[:,:,rec]

    @info "Extrema of Error for $(strDataset[l]): $(extrema(err))"
    @info "Mean of Error for $(strDataset[l]): $(mean(err))"

    ax = CairoMakie.Axis(fig[4, l], xlabel = rich("κ", subscript("y"), font=:italic), 
           ylabel = (l==1) ? rich("κ", subscript("x"), font=:italic) : "",
           title = "Error: $(strDataset[l])", xlabelpadding = 2, ylabelpadding = 2,
           yreversed = true, xticks = 1:2:maxMixingError, yticks = 1:2:maxMixingError,
           yticklabelsvisible=(l==1), yticksvisible=(l==1))

    @info maximum(err)
    CairoMakie.heatmap!(ax, transpose(err), colorrange = crange   )  
    tightlimits!(ax)
  end

  fig[4,5] = Colorbar(fig[4, 4], limits = crange, label = "Error")

  rowgap!(fig.layout, 10) 
  colgap!(fig.layout, 10)

  save("img/$(filenamePrefix)$(rec).png", fig, px_per_unit = 2, pt_per_unit = 0.25)
  save("img/$(filenamePrefix)$(rec).pdf", fig, px_per_unit = 2, pt_per_unit = 0.25)
end

