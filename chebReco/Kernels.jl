export DCR2_generate_kernels, DCR2_generate_aniso_kernels

function DCR2_generate_kernels(params, shape, basistype="UU" )

  DCore = params[:DCore]
  kB = 1.38064852e-23
  VCore = π/6 * DCore^3
  μ₀ = 4*π*1e-7 #vacuum permeability
  MS = 474000.0;
  temp = 293.0;
  msat = MS * VCore #saturation magnetic moment of a single nanoparticle
  beta = msat /(kB*temp) #H measured in T/μ₀


  area = params[:amplitude] ./ params[:gradient]   #0.012
  enlargementFactor = 2
  n_kern = (enlargementFactor) .* (shape.-1).+1;
  x2 = range(-enlargementFactor*area[1],enlargementFactor*area[1], length=n_kern[1] ); # larger evaluation of kernel to avoid artifacts
  y2 = range(-enlargementFactor*area[2],enlargementFactor*area[2], length=n_kern[2] );
  X2,Y2 = meshgrid(x2,y2);

  if basistype == "UU"

    function LangDxDy(x,y)
      if abs(x) < eps() && abs(y) < eps()
        return 0
      else 
        return (y*(-6.0*x^2 + 2.0 *y^2 + (2.0* x^4 + x^2 * y^2 - y^4) * 
              csch(sqrt(x^2 + y^2))^2 + sqrt(x^2 + y^2) * 
                coth(sqrt(x^2 + y^2)) * (2.0*x^2 - y^2 + 2.0* x^2 *(x^2 + y^2) * csch(sqrt(x^2 + y^2))^2)))/(x^2 + y^2)^3;
      end
    end
    global Kern_x = LangDxDy.(params[:gradient][1]*beta*vec(X2),params[:gradient][2]*beta*vec(Y2));
    Kern_x = reshape(Kern_x, n_kern[1:2]... );
    global Kern_y = LangDxDy.(params[:gradient][2]*beta*vec(Y2),params[:gradient][1]*beta*vec(X2));
    Kern_y = reshape(Kern_y, n_kern[1:2]... );

    return Kern_x, Kern_y

  elseif basistype == "UT"
  
    function LangDx_x(x,y) 
      if abs(x) < eps() && abs(y) < eps()
        return 1/3
      else
        return ( y^2 * sqrt(x^2 + y^2) * 
                coth(sqrt(x^2 + y^2)) + x^2 * (-(x^2 + y^2)) * csch(sqrt(x^2 + y^2))^2 + x^2 - y^2)/(x^2 + y^2)^2;
      end
    end
    function LangDy_y(x,y) 
      return LangDx_x(y,x)
    end

    global Kern_x = LangDx_x.(params[:gradient][1]*beta*vec(X2), params[:gradient][2]*beta*vec(Y2));
    Kern_x = reshape(Kern_x, n_kern[1:2]... );

    global Kern_y = LangDy_y.( params[:gradient][1]*beta*vec(X2),params[:gradient][2]*beta*vec(Y2));
    Kern_y = reshape(Kern_y, n_kern[1:2]... );

    return Kern_x, Kern_y
  end
end



function DCR2_generate_aniso_kernels(params, shape, split=1) #shape
  
  amplitude = params[:amplitude]
  gradient = params[:gradient] 
  grid = RegularGridPositions(collect((shape[1],shape[2],1)), 
            collect(2 .* (amplitude ./ gradient)),[0.0,0.0,0.0])  #params[:grid]

  offsets = vec([ gradient.*Tuple(x)  for x in grid ])

  numLenGrid = shape[1]*shape[2]

  factor = params[:kAnis]
  kAnisγ = params[:kAnisγ]
  maxOff = maximum(norm.(offsets))
  anisotropyAxis = vec([ (1.0./maxOff).*gradient.*(x[1],x[2],x[3])  for x in grid ])

  MS = 474000.0;
  temp = 293.0;
  
  
  # please do not change the sign of the selection field, otherwise the phase
  # must be manipulated

  GradStrx = gradient[1];
  GradStry = gradient[2];
  GradStrz = gradient[3];
  
  # amplitudes can be manipulated.
  Ax = amplitude[1]; # amplitude | x-direction drive field
  Ay = amplitude[2]; # amplitude | y-direction drive field
  Az = amplitude[3]; # amplitude | z-direction
 
  
  
  order = params[:order]; # order factor for truncation of series expansion -> large values could lead to NaNs or Inf.
  # small order faster calculation, but lower accuracy
  epsilon = 1e-10; #to avoid divisions by values close to zero in the series expansion
  
  
  nEasy_H_S = collect([ (norm(x) > 0 ? x[d]./norm(x) : (d==1)) for x in anisotropyAxis, d=1:3]')
  g_Aniso = factor.*(norm.(anisotropyAxis)).^kAnisγ
  
  enlargeFactor = 2
  
  shape_ = enlargeFactor .* (shape.-1).+1 
  
  #Create a Chebyshev grid around each sampling grid point.
  xk = enlargeFactor.*2.0*Ax./GradStrx.*(((-1:(shape_[1].-1)).+1/2 )./(shape_[1].-1) .-1/2);
  yk = enlargeFactor.*2.0*Ay./GradStry.*(((-1:(shape_[2].-1)).+1/2 )./(shape_[2].-1) .-1/2);
  
  # Coordinates that must be evaluated and padding of K_Aniso and the easy axis
  deltaX, deltaY = meshgrid(xk,yk);
  X_cords_gauss_quad = repeat(-deltaX, outer=[1,1,numLenGrid])
  Y_cords_gauss_quad = repeat(-deltaY, outer=[1,1,numLenGrid])
  Z_cords_gauss_quad = zeros(size(deltaX,1),size(deltaX,2),numLenGrid)

  # vectorized
  g_aniso_quad = vec(permutedims(repeat(g_Aniso,outer=[1,size(deltaX,1),size(deltaX,2)]),(2,3,1)))
  nEasy_quad = repeat(permutedims(reshape(nEasy_H_S,(size(nEasy_H_S)...,1,1)),(1,3,4,2)),outer=[1,size(deltaX,1),size(deltaX,2)])
  nEasy_quad = reshape(nEasy_quad,3,:)
  
  # The magentic field to be evaluated
  H = hcat(GradStrx*X_cords_gauss_quad[:],GradStry*Y_cords_gauss_quad[:],GradStrz*Z_cords_gauss_quad[:]);
  split += 1
  idx = Int.(ceil.(range(1,size(nEasy_quad,2),split)))
  idxS = idx[1:end-1]
  idxE = idx[2:end].-1
  idxE[end] = idx[end]
  
  Mag = zeros(size(nEasy_quad,2),3)
  for ind in range(1,length(idxS))
    idx = idxS[ind]:idxE[ind]
    H_red = H[idx,:]
    g_aniso_red = g_aniso_quad[idx]
    nEasy_red = nEasy_quad[:,idx]

    #calculated magnetic moment of the given field H
    Mag[idx,:] .= MNPDynamics.eqAnisoMeanMagFullVec(H_red, params[:DCore], MS, temp, g_aniso_red, 
                                            nEasy_red, order, epsilon);
  end
    
  Mag = reshape(Mag,size(X_cords_gauss_quad)..., 3);

  return diff(diff(Mag, dims=1),dims=2), diff(Mag, dims=2), diff(Mag, dims=1)
end
