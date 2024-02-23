
P_aufl1 = zeros(340,300)
P_aufl1[120:144,51:250] .= 1
P_aufl1[179:203,63:238] .= 1
P_aufl1_O = zeros(17,15)

for x=1:340,y=1:300
    P_aufl1_O[Int(ceil(x/20)),Int(ceil(y/20))] += P_aufl1[x,y]/400#itp(-10+x*20,-10+y*20)
end

## Pos 2
# Balken 1 in 7. Pixelreihe -> 120-144, mittig -> 51-250
# Balken 2 in 10.-11. Pixelreihe & ~5mm (50 Reihen) unter Balken 1, mittig -> 63-238

P_aufl2 = zeros(340,300)
P_aufl2[120:144,51:250] .= 1
P_aufl2[190:214,63:238] .= 1
P_aufl2_O = zeros(17,15)

for x=1:340,y=1:300
    P_aufl2_O[Int(ceil(x/20)),Int(ceil(y/20))] += P_aufl2[x,y]/400#itp(-10+x*20,-10+y*20)
end


## Pos 3
# Balken 1 in 7. Pixelreihe -> 120-144, mittig -> 51-250
# Balken 3 in 11.-12. Pixelreihe & ~7mm (70 Reihen) unter Balken 1, mittig -> 63-238

P_aufl3 = zeros(340,300)
P_aufl3[120:144,51:250] .= 1
P_aufl3[209:233,63:238] .= 1
P_aufl3_O = zeros(17,15)

for x=1:340,y=1:300
    P_aufl3_O[Int(ceil(x/20)),Int(ceil(y/20))] += P_aufl3[x,y]/400#itp(-10+x*20,-10+y*20)
end

## SpiralPhantom

P_spiral = zeros(340,300)
P_spiral[65:264,46:70] .= 1
P_spiral[42:66,84:258] .= 1
P_spiral[80:229,239:263] .= 1
P_spiral[232:256,141:227] .= 1
P_spiral[171:220,136:161] .= 1
P_spiral_O = zeros(17,15)

for x=1:340,y=1:300
    P_spiral_O[Int(ceil(x/20)),Int(ceil(y/20))] += P_spiral[x,y]/400#itp(-10+x*20,-10+y*20)
end

## ICEPhantom

function myhalvecircle(N,M;rad=1.0,mid=(1+N,1+M)./2)
	Inds = CartesianIndex[]
    for ind in mycircle(N,M;rad=rad,mid=mid)
        if ind[1] <= mid[1]
            push!(Inds,ind)
        end
    end
    return Inds
end


function mycircle(N,M;rad=1.0,mid=(1+N,1+M)./2)
	return myellipse(N,M;a=rad,b=rad,mid=mid)
end

function myellipse(N,M;a=1.0,b=1.0,mid=(1+N,1+M)./2)
	allInds = CartesianIndex(1,1):CartesianIndex(N,M)
	Inds = CartesianIndex[]
	for ind in allInds
		pos = ind.I .- mid
		if pos[1]^2/a^2 + pos[2]^2/b^2 <= 1
			push!(Inds,ind)
		end
	end
	return Inds
end

function mycone(N,M;h=1,r=1,mid=(1+N,1+M)./2)
	allInds = CartesianIndex(1,1):CartesianIndex(N,M)
	Inds = CartesianIndex[]
	for ind in allInds
		pos = ind.I .- mid
		if 0 <= pos[1] <= h
            if abs(pos[2]) <= (1-pos[1]/h) * r
			push!(Inds,ind)
            end
		end
	end
	return Inds
end

P_ICE = zeros(340,300)
for ind in myhalvecircle(340,300;rad=65,mid=(135,150.5))
    P_ICE[ind] =1
end
for ind in mycone(340,300;h=112,r=65,mid=(136,150.5))
    P_ICE[ind] =1
end
P_ICE_O = zeros(17,15)

for x=1:340,y=1:300
    P_ICE_O[Int(ceil(x/20)),Int(ceil(y/20))] += P_ICE[x,y]/400#itp(-10+x*20,-10+y*20)
end

# dot 6mm in x and y from mid (171,151)+(60,60)=(231,211)+-(10,10)
P_dot = zeros(340,300)
P_dot[221:240,201:220] .= 1
P_dot_O = zeros(17,15)

for x=1:340,y=1:300
    P_dot_O[Int(ceil(x/20)),Int(ceil(y/20))] += P_dot[x,y]/400#itp(-10+x*20,-10+y*20)
end



DIE_MASKE=reverse.(transpose.([P_spiral_O,P_aufl1_O,P_aufl2_O,P_aufl3_O,P_ICE_O,P_dot_O]),dims=2)
