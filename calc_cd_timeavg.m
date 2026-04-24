function cdImg = calc_cd_timeavg(MAG, vMean, Venc)
%CALC_CD_TIMEAVG  PC-style CD angiogram from time-averaged MAG and velocity (not temporal MIP).
%
%   cdImg = MAG .* sin( (pi/2) * min(|v|, Venc) / Venc )
%
%   MAG: 3-D magnitude; vMean: 4-D [nx,ny,nz,3] mean velocity; Venc scalar (same units as v).

MAG = single(MAG);
vMean = single(vMean);
Vmag = sqrt(sum(vMean.^2, 4));
Vmag = min(Vmag, single(Venc));
Vmag = max(Vmag, 0);
vencS = max(single(Venc), single(1e-6));
cdImg = MAG .* sin((single(pi) / 2) .* (Vmag ./ vencS));
cdImg = max(cdImg, 0);

end
