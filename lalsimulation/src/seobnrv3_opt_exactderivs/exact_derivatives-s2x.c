const double sigmaKerrdata0prm  =  1.0;
const double sigmaStardata0prm  =  mass1overmass2;
const double s2dots2prm  =  2*s2Vec->data[0];
const double sKerrUSCORExprm  =  sigmaKerrdata0prm;
const double sStarUSCORExprm  =  sigmaStardata0prm;
const double a2prm  =  2*sKerrUSCOREx*sKerrUSCORExprm;
const double a4prm  =  2*a2*a2prm;
const double aprm  =  a2prm/(2.*Sqrt(a2));
const double coeffsk2prm  =  c1k2*a2prm;
const double coeffsk3prm  =  c1k3*a2prm;
const double coeffsk4prm  =  c1k4*a2prm + c2k4*a4prm;
const double coeffsk5prm  =  c1k5*a2prm + c2k5*a4prm;
const double costhetaprm  =  nx*e3USCORExprm + ny*e3USCOREyprm + nz*e3USCOREzprm;
const double xi2prm  =  -2*costheta*costhetaprm;
const double xiUSCORExprm  =  nz*e3USCOREyprm - ny*e3USCOREzprm;
const double xiUSCOREyprm  =  -(nz*e3USCORExprm) + nx*e3USCOREzprm;
const double xiUSCOREzprm  =  ny*e3USCORExprm - nx*e3USCOREyprm;
const double vxprm  =  -(nz*xiUSCOREyprm) + ny*xiUSCOREzprm;
const double vyprm  =  nz*xiUSCORExprm - nx*xiUSCOREzprm;
const double vzprm  =  -(ny*xiUSCORExprm) + nx*xiUSCOREyprm;
const double w2prm  =  a2prm;
const double rho2prm  =  costheta*(costheta*a2prm + 2*a2*costhetaprm);
const double bulkprm  =  u2*a2prm;
const double logargprm  =  u2*coeffsk2prm + u3*coeffsk3prm + u4*coeffsk4prm + u5*coeffsk5prm;
const double onepluslogargprm  =  logargprm;
const double invonepluslogargprm  =  (-onepluslogargprm)/((onepluslogarg)*(onepluslogarg));
const double logTermsprm  =  (eta*onepluslogargprm)/onepluslogarg;
const double deltaUprm  =  logTerms*bulkprm + bulk*logTermsprm;
const double deltaTprm  =  r2*deltaUprm;
const double deltaUUSCOREupt7prm  =  coeffsk5prm;
const double deltaUUSCOREupt6prm  =  4.*coeffsk4prm + 5.*u*deltaUUSCOREupt7prm;
const double deltaUUSCOREupt5prm  =  3.*coeffsk3prm + u*deltaUUSCOREupt6prm;
const double deltaUUSCOREupt4prm  =  2.*coeffsk2prm + u*deltaUUSCOREupt5prm;
const double deltaUUSCOREupt3prm  =  u*deltaUUSCOREupt4prm;
const double deltaUUSCOREupt2prm  =  u*a2prm;
const double deltaUUSCOREupt1prm  =  eta*(deltaUUSCOREupt3*bulkprm + bulk*deltaUUSCOREupt3prm);
const double deltaUUSCOREuprm  =  invonepluslogarg*deltaUUSCOREupt1prm + 2.*logTerms*deltaUUSCOREupt2prm + deltaUUSCOREupt1*invonepluslogargprm + 2.*deltaUUSCOREupt2*logTermsprm;
const double deltaTUSCORErprm  =  2.*r*deltaUprm - deltaUUSCOREuprm;
const double Lambdaprm  =  -(a2*xi2*deltaTprm) + 2*w2*w2prm - deltaT*(xi2*a2prm + a2*xi2prm);
const double rho2xi2Lambdaprm  =  Lambda*xi2*rho2prm + rho2*(xi2*Lambdaprm + Lambda*xi2prm);
const double invrho2xi2Lambdaprm  =  (-rho2xi2Lambdaprm)/((rho2xi2Lambda)*(rho2xi2Lambda));
const double invrho2prm  =  invrho2xi2Lambda*xi2*Lambdaprm + Lambda*(xi2*invrho2xi2Lambdaprm + invrho2xi2Lambda*xi2prm);
const double invxi2prm  =  invrho2xi2Lambda*rho2*Lambdaprm + Lambda*(rho2*invrho2xi2Lambdaprm + invrho2xi2Lambda*rho2prm);
const double invLambdaprm  =  invrho2xi2Lambda*xi2*rho2prm + rho2*(xi2*invrho2xi2Lambdaprm + invrho2xi2Lambda*xi2prm);
const double invLambdasqprm  =  2*invLambda*invLambdaprm;
const double rho2invLambdaprm  =  rho2*invLambdaprm + invLambda*rho2prm;
const double expnuprm  =  (rho2invLambda*deltaTprm + deltaT*rho2invLambdaprm)/(2.*Sqrt(deltaT*rho2invLambda));
const double expMUprm  =  rho2prm/(2.*Sqrt(rho2));
const double expMUexpnuprm  =  expnu*expMUprm + expMU*expnuprm;
const double expMUsqprm  =  2*expMU*expMUprm;
const double expnusqprm  =  2*expnu*expnuprm;
const double expMUsqexpnusqprm  =  expnusq*expMUsqprm + expMUsq*expnusqprm;
const double invexpnuexpMUprm  =  (-expMUexpnuprm)/((expMUexpnu)*(expMUexpnu));
const double invexpMUprm  =  invexpnuexpMU*expnuprm + expnu*invexpnuexpMUprm;
const double invexpMUsqprm  =  2*invexpMU*invexpMUprm;
const double expnuinvexpMU2prm  =  invexpMUsq*expnuprm + expnu*invexpMUsqprm;
const double invexpMUcubinvexpnuprm  =  invexpnuexpMU*invexpMUsqprm + invexpMUsq*invexpnuexpMUprm;
const double deltaRprm  =  DD*deltaTprm;
const double wwprm  =  (2.*r + coeffs->bb3*eta*u + coeffs->b3*eta*u*a2)*aprm + coeffs->b3*eta*u*a*a2prm;
const double Bprm  =  deltaTprm/(2.*Sqrt(deltaT));
const double sqrtdeltaTprm  =  Bprm;
const double sqrtdeltaRprm  =  deltaRprm/(2.*Sqrt(deltaR));
const double deltaTsqrtdeltaRprm  =  sqrtdeltaR*deltaTprm + deltaT*sqrtdeltaRprm;
const double sqrtdeltaTdeltaTsqrtdeltaRprm  =  sqrtdeltaT*deltaTsqrtdeltaRprm + deltaTsqrtdeltaR*sqrtdeltaTprm;
const double invdeltaTsqrtdeltaTsqrtdeltaRprm  =  (-sqrtdeltaTdeltaTsqrtdeltaRprm)/((sqrtdeltaTdeltaTsqrtdeltaR)*(sqrtdeltaTdeltaTsqrtdeltaR));
const double invdeltaTprm  =  invdeltaTsqrtdeltaTsqrtdeltaR*sqrtdeltaT*sqrtdeltaRprm + sqrtdeltaR*(sqrtdeltaT*invdeltaTsqrtdeltaTsqrtdeltaRprm + invdeltaTsqrtdeltaTsqrtdeltaR*sqrtdeltaTprm);
const double invsqrtdeltaTprm  =  invdeltaTsqrtdeltaTsqrtdeltaR*deltaTsqrtdeltaRprm + deltaTsqrtdeltaR*invdeltaTsqrtdeltaTsqrtdeltaRprm;
const double invsqrtdeltaRprm  =  deltaT*sqrtdeltaT*invdeltaTsqrtdeltaTsqrtdeltaRprm + invdeltaTsqrtdeltaTsqrtdeltaR*(sqrtdeltaT*deltaTprm + deltaT*sqrtdeltaTprm);
const double wprm  =  ww*invLambdaprm + invLambda*wwprm;
const double LambdaUSCORErprm  =  -(a2*xi2*deltaTUSCORErprm) + 4.*r*w2prm - deltaTUSCOREr*(xi2*a2prm + a2*xi2prm);
const double wwUSCORErprm  =  -(-2. + coeffs->bb3*eta*u2 + coeffs->b3*eta*u2*a2)*aprm - coeffs->b3*eta*u2*a*a2prm;
const double BRprm  =  -invsqrtdeltaT*(invsqrtdeltaR*deltaTprm - 0.5*deltaTUSCORErprm + deltaT*invsqrtdeltaRprm) + (0.5*deltaTUSCOREr - deltaT*invsqrtdeltaR)*invsqrtdeltaTprm;
const double wrprm  =  (-(LambdaUSCOREr*ww) + Lambda*wwUSCOREr)*invLambdasqprm + invLambdasq*(wwUSCOREr*Lambdaprm - ww*LambdaUSCORErprm - LambdaUSCOREr*wwprm + Lambda*wwUSCORErprm);
const double nurpt2prm  =  ((w2)*(w2))*deltaTUSCORErprm - 4.*r*deltaT*w2prm + w2*(-4.*r*deltaTprm + 2*deltaTUSCOREr*w2prm);
const double nurpt1prm  =  nurpt2*invdeltaTprm + invdeltaT*nurpt2prm;
const double nurprm  =  0.5*nurpt1*invLambdaprm + r*invrho2prm + 0.5*invLambda*nurpt1prm;
const double murprm  =  r*invrho2prm - invsqrtdeltaRprm;
const double a2costhetaprm  =  costheta*a2prm + a2*costhetaprm;
const double wcospt2prm  =  ww*deltaTprm + deltaT*wwprm;
const double wcospt1prm  =  wcospt2*invLambdasqprm + invLambdasq*wcospt2prm;
const double wcosprm  =  -2.*wcospt1*a2costhetaprm - 2.*a2costheta*wcospt1prm;
const double nucospt3prm  =  invrho2*invLambdaprm + invLambda*invrho2prm;
const double nucospt2prm  =  w2*nucospt3prm + nucospt3*w2prm;
const double nucospt1prm  =  nucospt2*a2costhetaprm + a2costheta*nucospt2prm;
const double nucosprm  =  (-deltaT + w2)*nucospt1prm + nucospt1*(-deltaTprm + w2prm);
const double mucosprm  =  invrho2*a2costhetaprm + a2costheta*invrho2prm;
const double csiprm  =  (deltaR*w2*deltaTprm + deltaT*(w2*deltaRprm - 2*deltaR*w2prm))/(2.*Sqrt(deltaR*deltaT)*((w2)*(w2)));
const double csi1prm  =  (1.-fabs(1.-tortoise)) * (csiprm);
const double csi2prm  =  (0.5-copysign(0.5, 1.5-tortoise)) * (csiprm);
const double prTprm  =  (nx*p->data[0] + ny*p->data[1] + nz*p->data[2])*csi2prm;
const double prTtimesoneminuscsi1invprm  =  (prT*csi1prm)/((csi1)*(csi1)) + (1. - 1./csi1)*prTprm;
const double tmpP0prm  =  -(nx*prTtimesoneminuscsi1invprm);
const double tmpP1prm  =  -(ny*prTtimesoneminuscsi1invprm);
const double tmpP2prm  =  -(nz*prTtimesoneminuscsi1invprm);
const double pxirprm  =  r*(xiUSCOREx*tmpP0prm + xiUSCOREy*tmpP1prm + xiUSCOREz*tmpP2prm + tmpP0*xiUSCORExprm + tmpP1*xiUSCOREyprm + tmpP2*xiUSCOREzprm);
const double pvrprm  =  r*(vx*tmpP0prm + vy*tmpP1prm + vz*tmpP2prm + tmpP0*vxprm + tmpP1*vyprm + tmpP2*vzprm);
const double pvrsqprm  =  2*pvr*pvrprm;
const double pnprm  =  nx*tmpP0prm + ny*tmpP1prm + nz*tmpP2prm;
const double pnsqprm  =  2*pn*pnprm;
const double prprm  =  pnprm;
const double prsqprm  =  2*pr*prprm;
const double pfprm  =  pxirprm;
const double pxirsqprm  =  2*pxir*pxirprm;
const double ptheta2prm  =  pvrsq*invxi2prm + invxi2*pvrsqprm;
const double prT4prm  =  4*((prT)*(prT)*(prT))*prTprm;
const double Hnspt7prm  =  invrho2*deltaRprm + deltaR*invrho2prm;
const double Hnspt6prm  =  rho2invLambda*invxi2prm + invxi2*rho2invLambdaprm;
const double Hnspt4prm  =  ((pf)*(pf))*Hnspt6prm + prsq*Hnspt7prm + ptheta2*invrho2prm + 2*Hnspt6*pf*pfprm + Hnspt7*prsqprm + Hnspt5*prT4prm + invrho2*ptheta2prm;
const double Hnspt3prm  =  Hnspt4*deltaTprm + deltaT*Hnspt4prm;
const double Hnspt2prm  =  rho2*Hnspt3prm + Hnspt3*rho2prm;
const double Hnspt1prm  =  ww*pfprm + pf*wwprm;
const double Hnsprm  =  invLambda*Hnspt1prm + Hnspt1*invLambdaprm + (invLambda*Hnspt2prm + Hnspt2*invLambdaprm)/(2.*Sqrt(Hnspt2*invLambda));
const double Qpt3prm  =  invrho2*deltaRprm + deltaR*invrho2prm;
const double Qpt2prm  =  rho2invLambda*invxi2prm + invxi2*rho2invLambdaprm;
const double Qpt1prm  =  invxi2*invrho2prm + invrho2*invxi2prm;
const double Qprm  =  Qpt3*pnsqprm + Qpt1*pvrsqprm + Qpt2*pxirsqprm + pvrsq*Qpt1prm + pxirsq*Qpt2prm + pnsq*Qpt3prm;
const double pn2prm  =  deltaR*prsq*invrho2prm + invrho2*(prsq*deltaRprm + deltaR*prsqprm);
const double ppprm  =  Qprm;
const double sKerrmultfactprm  =  -36.*r*pn2prm + 3.*r*ppprm;
const double sStarmultfactprm  =  r*(-30.*pn2prm + 4.*ppprm);
const double deltaSigmaStarUSCOREx1prm  =  etaover12r*(sKerrUSCOREx*sKerrmultfactprm + sKerrmultfact*sKerrUSCORExprm + sStarUSCOREx*sStarmultfactprm + sStarmultfact*sStarUSCORExprm);
const double deltaSigmaStarUSCOREy1prm  =  etaover12r*(sKerrUSCOREy*sKerrmultfactprm + sStarUSCOREy*sStarmultfactprm);
const double deltaSigmaStarUSCOREz1prm  =  etaover12r*(sKerrUSCOREz*sKerrmultfactprm + sStarUSCOREz*sStarmultfactprm);
const double pn2ppprm  =  pp*pn2prm + pn2*ppprm;
const double pp2prm  =  2*pp*ppprm;
const double pn2u2prm  =  u2*pn2prm;
const double ppu2prm  =  u2*ppprm;
const double pn2ppu2prm  =  u2*pn2ppprm;
const double sMultiplier1pt6prm  =  -720.*pn2*pn2prm + 126.*pn2ppprm + 3.*pp2prm;
const double sMultiplier1pt5prm  =  -96.*pn2ppprm + 23.*pp2prm;
const double sMultiplier1pt4prm  =  324.*pn2prm - 120.*ppprm + r*sMultiplier1pt6prm;
const double sMultiplier1pt3prm  =  -282.*pn2prm + 206.*ppprm + r*sMultiplier1pt5prm;
const double sMultiplier1pt2prm  =  r*sMultiplier1pt4prm;
const double sMultiplier1pt1prm  =  eta*sMultiplier1pt2prm + r*sMultiplier1pt3prm;
const double sMultiplier1prm  =  -0.013888888888888888*eta*u2*sMultiplier1pt1prm;
const double sMultiplier2pt6prm  =  5.625*pn2u2*pn2prm - 1.625*pn2ppu2prm + 5.625*pn2*pn2u2prm;
const double sMultiplier2pt5prm  =  0.25*pn2ppu2prm - 0.3125*u2*pp2prm;
const double sMultiplier2pt4prm  =  -6.125*pn2u2prm + 1.4166666666666665*ppu2prm + r*sMultiplier2pt6prm;
const double sMultiplier2pt3prm  =  -0.6666666666666666*pn2u2prm - 3.0277777777777777*ppu2prm + r*sMultiplier2pt5prm;
const double sMultiplier2pt2prm  =  r*sMultiplier2pt4prm;
const double sMultiplier2pt1prm  =  eta*sMultiplier2pt2prm + r*sMultiplier2pt3prm;
const double sMultiplier2prm  =  eta*sMultiplier2pt1prm;
const double deltaSigmaStarUSCOREx2prm  =  deltaSigmaStarUSCOREx1prm + sMultiplier2*sigmaKerrdata0prm + sMultiplier1*sigmaStardata0prm + sigmaStar->data[0]*sMultiplier1prm + sigmaKerr->data[0]*sMultiplier2prm;
const double deltaSigmaStarUSCOREy2prm  =  deltaSigmaStarUSCOREy1prm + sigmaStar->data[1]*sMultiplier1prm + sigmaKerr->data[1]*sMultiplier2prm;
const double deltaSigmaStarUSCOREz2prm  =  deltaSigmaStarUSCOREz1prm + sigmaStar->data[2]*sMultiplier1prm + sigmaKerr->data[2]*sMultiplier2prm;
const double deltaSigmaStarUSCOREx3prm  =  deltaSigmaStarUSCOREx2prm + coeffs->d1*etau3*sigmaStardata0prm;
const double deltaSigmaStarUSCOREy3prm  =  deltaSigmaStarUSCOREy2prm;
const double deltaSigmaStarUSCOREz3prm  =  deltaSigmaStarUSCOREz2prm;
const double deltaSigmaStarUSCORExprm  =  deltaSigmaStarUSCOREx3prm + coeffs->d1v2*etau3*sigmaKerrdata0prm;
const double deltaSigmaStarUSCOREyprm  =  deltaSigmaStarUSCOREy3prm;
const double deltaSigmaStarUSCOREzprm  =  deltaSigmaStarUSCOREz3prm;
const double sxprm  =  deltaSigmaStarUSCORExprm + sStarUSCORExprm;
const double syprm  =  deltaSigmaStarUSCOREyprm;
const double szprm  =  deltaSigmaStarUSCOREzprm;
const double sxiprm  =  xiUSCOREx*sxprm + xiUSCOREy*syprm + xiUSCOREz*szprm + sx*xiUSCORExprm + sy*xiUSCOREyprm + sz*xiUSCOREzprm;
const double svprm  =  vx*sxprm + vy*syprm + vz*szprm + sx*vxprm + sy*vyprm + sz*vzprm;
const double snprm  =  nx*sxprm + ny*syprm + nz*szprm;
const double s3prm  =  sx*e3USCORExprm + sy*e3USCOREyprm + sz*e3USCOREzprm + e3USCOREx*sxprm + e3USCOREy*syprm + e3USCOREz*szprm;
const double sqrtQprm  =  Qprm/(2.*Sqrt(Q));
const double oneplus2sqrtQprm  =  2.*sqrtQprm;
const double oneplus1sqrtQprm  =  oneplus2sqrtQprm - sqrtQprm;
const double twoB1psqrtQsqrtQprm  =  2.*B*sqrtQ*oneplus1sqrtQprm + oneplus1sqrtQ*(2.*sqrtQ*Bprm + 2.*B*sqrtQprm);
const double invtwoB1psqrtQsqrtQprm  =  (-twoB1psqrtQsqrtQprm)/((twoB1psqrtQsqrtQ)*(twoB1psqrtQsqrtQ));
const double expMUsqsqrtQplusQprm  =  (Q + sqrtQ)*expMUsqprm + expMUsq*(Qprm + sqrtQprm);
const double Hwrpt4aprm  =  sv*pxirsqprm + pxirsq*svprm;
const double Hwrpt4prm  =  Hwrpt4a*expMUsqexpnusqprm + expMUsqexpnusq*Hwrpt4aprm;
const double Hwrpt3cprm  =  sxi*pxirprm + pxir*sxiprm;
const double Hwrpt3bprm  =  pvr*Hwrpt3cprm + Hwrpt3c*pvrprm;
const double Hwrpt3aprm  =  Hwrpt3b*expMUexpnuprm + expMUexpnu*Hwrpt3bprm;
const double Hwrpt3prm  =  Hwrpt3a*Bprm + B*Hwrpt3aprm;
const double Hwrpt2gprm  =  sv*deltaRprm + deltaR*svprm;
const double Hwrpt2fprm  =  sqrtdeltaR*snprm + sn*sqrtdeltaRprm;
const double Hwrpt2eprm  =  pvr*Hwrpt2fprm + Hwrpt2f*pvrprm;
const double Hwrpt2dprm  =  pnsq*Hwrpt2gprm + Hwrpt2g*pnsqprm;
const double Hwrpt2cprm  =  pn*Hwrpt2eprm + Hwrpt2e*pnprm;
const double Hwrpt2bprm  =  sv*expMUsqsqrtQplusQprm + expMUsqsqrtQplusQ*svprm;
const double Hwrpt2aprm  =  xi2*(Hwrpt2bprm + Hwrpt2cprm - Hwrpt2dprm) + (Hwrpt2b + Hwrpt2c - Hwrpt2d)*xi2prm;
const double Hwrpt2prm  =  Hwrpt2a*deltaTprm + deltaT*Hwrpt2aprm;
const double Hwrpt1bprm  =  invxi2*invtwoB1psqrtQsqrtQprm + invtwoB1psqrtQsqrtQ*invxi2prm;
const double Hwrpt1aprm  =  sqrtdeltaR*Hwrpt1bprm + Hwrpt1b*sqrtdeltaRprm;
const double Hwrpt1prm  =  invexpMUcubinvexpnu*Hwrpt1aprm + Hwrpt1a*invexpMUcubinvexpnuprm;
const double Hwrprm  =  (Hwrpt2 - Hwrpt3 + Hwrpt4)*Hwrpt1prm + Hwrpt1*(Hwrpt2prm - Hwrpt3prm + Hwrpt4prm);
const double Hwcospt9prm  =  sxi*pxirprm + pxir*sxiprm;
const double Hwcospt8prm  =  sv*pvrprm + pvr*svprm;
const double Hwcospt7prm  =  Hwcospt8*Bprm - Hwcospt9*expMUexpnuprm + B*Hwcospt8prm - expMUexpnu*Hwcospt9prm;
const double Hwcospt6prm  =  sqrtdeltaR*Hwcospt7prm + Hwcospt7*sqrtdeltaRprm;
const double Hwcospt5prm  =  -(xi2*expMUsqsqrtQplusQprm) + pvrsqprm - expMUsqsqrtQplusQ*xi2prm;
const double Hwcospt4prm  =  pn*Hwcospt6prm + Hwcospt6*pnprm;
const double Hwcospt3prm  =  Hwcospt5*deltaTprm - pxirsq*expMUsqexpnusqprm + deltaT*Hwcospt5prm - expMUsqexpnusq*pxirsqprm;
const double Hwcospt2prm  =  -(Hwcospt4*Bprm) + sn*Hwcospt3prm - B*Hwcospt4prm + Hwcospt3*snprm;
const double Hwcospt1prm  =  invexpMUcubinvexpnu*Hwcospt2prm + Hwcospt2*invexpMUcubinvexpnuprm;
const double Hwcosprm  =  invtwoB1psqrtQsqrtQ*Hwcospt1prm + Hwcospt1*invtwoB1psqrtQsqrtQprm;
const double deltaTsqrtQprm  =  sqrtQ*deltaTprm + deltaT*sqrtQprm;
const double invdeltatTsqrtQprm  =  (-deltaTsqrtQprm)/((deltaTsqrtQ)*(deltaTsqrtQ));
const double HSOLpt5prm  =  pxir*(-Bprm + expMUexpnuprm) + (-B + expMUexpnu)*pxirprm;
const double HSOLpt4prm  =  invexpMU*HSOLpt5prm + HSOLpt5*invexpMUprm;
const double HSOLpt3prm  =  HSOLpt4*expnusqprm + expnusq*HSOLpt4prm;
const double HSOLpt2prm  =  s3*HSOLpt3prm + HSOLpt3*s3prm;
const double HSOLpt1prm  =  invxi2*HSOLpt2prm + HSOLpt2*invxi2prm;
const double HSOLprm  =  invdeltatTsqrtQ*HSOLpt1prm + HSOLpt1*invdeltatTsqrtQprm;
const double deltaTsqrtQplusQprm  =  (Q + sqrtQ)*deltaTprm + deltaT*(Qprm + sqrtQprm);
const double invdeltaTsqrtQplusQprm  =  (-deltaTsqrtQplusQprm)/((deltaTsqrtQplusQ)*(deltaTsqrtQplusQ));
const double HSONLmult2prm  =  invxi2*invdeltaTsqrtQplusQprm + invdeltaTsqrtQplusQ*invxi2prm;
const double HSONLmultprm  =  HSONLmult2*expnuinvexpMU2prm + expnuinvexpMU2*HSONLmult2prm;
const double HSONLpt1bprm  =  xi2*pnprm + pn*xi2prm;
const double HSONLpt1aprm  =  (-mucos + nucos)*HSONLpt1bprm + pvr*murprm + HSONLpt1b*(-mucosprm + nucosprm) - pvr*nurprm + mur*pvrprm - nur*pvrprm;
const double HSONLpt1prm  =  sqrtQ*HSONLpt1aprm - mucos*HSONLpt1bprm - HSONLpt1b*mucosprm + pvr*murprm + mur*pvrprm + HSONLpt1a*sqrtQprm;
const double HSONLpt2dprm  =  pxir*nurprm + nur*pxirprm;
const double HSONLpt2cprm  =  oneplus2sqrtQ*HSONLpt2dprm + HSONLpt2d*oneplus2sqrtQprm;
const double HSONLpt2bprm  =  sxi*Bprm + B*sxiprm;
const double HSONLpt2aprm  =  HSONLpt2c*expMUexpnuprm + expMUexpnu*HSONLpt2cprm;
const double HSONLpt2prm  =  HSONLpt2b*HSONLpt1prm + sv*HSONLpt2aprm + HSONLpt1*HSONLpt2bprm + HSONLpt2a*svprm;
const double HSONLpt3cprm  =  sv*pxirprm + pxir*svprm;
const double HSONLpt3bprm  =  oneplus1sqrtQ*HSONLpt3cprm + HSONLpt3c*oneplus1sqrtQprm;
const double HSONLpt3aprm  =  HSONLpt3b*expMUexpnuprm + expMUexpnu*HSONLpt3bprm;
const double HSONLpt3prm  =  HSONLpt2*Bprm - HSONLpt3a*BRprm + B*HSONLpt2prm - BR*HSONLpt3aprm;
const double HSONLpt4eprm  =  xi2*snprm + sn*xi2prm;
const double HSONLpt4dprm  =  oneplus2sqrtQ*HSONLpt4eprm + HSONLpt4e*oneplus2sqrtQprm;
const double HSONLpt4cprm  =  pxir*HSONLpt4dprm + HSONLpt4d*pxirprm;
const double HSONLpt4bprm  =  nucos*HSONLpt4cprm + HSONLpt4c*nucosprm;
const double HSONLpt4aprm  =  HSONLpt4b*expMUexpnuprm + expMUexpnu*HSONLpt4bprm;
const double HSONLpt4prm  =  -(HSONLpt4a*Bprm) + sqrtdeltaR*HSONLpt3prm - B*HSONLpt4aprm + HSONLpt3*sqrtdeltaRprm;
const double HSONLprm  =  HSONLpt4*HSONLmultprm + HSONLmult*HSONLpt4prm;
const double Hsprm  =  HSOLprm + HSONLprm + wcos*Hwcosprm + wr*Hwrprm + w*s3prm + s3*wprm + Hwcos*wcosprm + Hwr*wrprm;
const double Hsspt1prm  =  -0.5*(-6.*sn*snprm + 2*(sx*sxprm + sy*syprm + sz*szprm));
const double Hssprm  =  u3*Hsspt1prm;
const double sKerrdotsStarprm  =  sStarUSCOREx*sKerrUSCORExprm + sKerrUSCOREx*sStarUSCORExprm;
const double Hpt1prm  =  etau4*s2dots2prm;
const double Hprm  =  Hnsprm + (coeffs->dheffSSv2 + coeffs->dheffSS*sKerrdotsStar)*Hpt1prm + Hsprm + Hssprm + coeffs->dheffSS*Hpt1*sKerrdotsStarprm;
const double Hrealprm  =  (eta*Hprm)*invHreal;