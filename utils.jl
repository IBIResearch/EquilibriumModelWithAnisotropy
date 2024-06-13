# Transfer function correction / estimation


# Estimate the transfer function of the receive chain using linear regression
function estimateTransferFunction(SMeas, SModel)
  a = zeros(ComplexF64, size(SMeas,3), size(SMeas,4))
  for k in CartesianIndices(a)
    a[k] = sum(conj(vec(SModel[:,:,k])) .* vec(SMeas[:,:,k])) / norm(vec(SModel[:,:,k])).^2
  end
  return a
end

function correctTransferFunction(SMeas, SModel)
  a = estimateTransferFunction(SMeas, SModel)
  SModel = reshape(a,1,1,size(SModel,3),size(SModel,4)) .* SModel
  SModel[:,:,1,:] .= 1
  return SModel
end


function estimateTransferFunction(SMeas, SModel, tfMeasured, shift, phi)
  transfer_signs = [1im,1im,1im]

  # shift correction
  N_ = (size(SModel,3)-1)*2
  shiftTerm = reshape( exp.(2*pi*im*phi .-2*pi*im.*(0:(size(SModel,3)-1))/N_*shift ),:,1)

  # here we to the scaling using mixing factor mx=my=3
  fac = maximum(abs.(SMeas[:, :, (3-1)*16+(3-1)*17+1, :]),dims=(1,2)) ./ 
      maximum(abs.(SModel[:, :, (3-1)*16+(3-1)*17+1, :]),dims=(1,2))

  # build the TF
  a_ = fac[1,:,:] .* (((tfMeasured .* reshape(transfer_signs,1,:))  ) .* shiftTerm) 

  return a_
end

function correctTransferFunction(SMeas, SModel, tfMeasured, shift, phi)
  a = estimateTransferFunction(SMeas, SModel,tfMeasured, shift, phi)
  SModel = reshape(a,1,1,size(SModel,3),size(SModel,4)) .* SModel
  SModel[:,:,1,:] .= 1
  return SModel
end

# Error calculations

function relErrorMixingOrder(smEstimate, smOrig, mixingOrders)
  err = zeros(length(mixingOrders),length(mixingOrders),2)

  for (x,my) in enumerate(mixingOrders)
    for (y,mx) in enumerate(mixingOrders)
      for r = 1:2

        A = smEstimate[:,:,(mx-1)*16+(my-1)*17+1, r]
        B = smOrig[:,:,(mx-1)*16+(my-1)*17+1, r]

        err[y, x, r] = norm(A-B) / maximum(abs.(B)) / sqrt(length(B))
      end
    end
  end

  if mixingOrders[1] == 1
    err[1, 1, :] .= err[2, 1, :] .= err[1, 2, :] .= 0 
  end
  return err
end

# colors

colors = [(0/255,73/255,146/255), # blue
          (239/255,123/255,5/255),	# orange (dark)
          (138/255,189/255,36/255),	# green
          (178/255,34/255,41/255), # red
          (170/255,156/255,143/255), # mocca
          (87/255,87/255,86/255),	# black (text)
          (255/255,223/255,0/255), # yellow
          (104/255,195/255,205/255),# "TUHH"
          (45/255,198/255,214/255), #  TUHH
          (193/255,216/255,237/255)]


function meshgrid(x, y)
  X = [x for _ in y, x in x]
  Y = [y for y in y, _ in x]
  X, Y
end

function meshgrid(x, y, z)
  X = [x for _ in z, _ in y, x in x]
  Y = [y for _ in z, y in y, _ in x]
  Z = [z for _ in z, y in y, _ in x]
  X, Y, Z
end