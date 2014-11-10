function noiseVar = aco_clipping_noise(signalVar,bottom,top)
lamdaBottom = bottom/sqrt(signalVar);
lamdaTop = top/sqrt(signalVar);
K = qfunc(lamdaBottom) - qfunc(lamdaTop);
noiseVar = K*(lamdaBottom^2+1)-2*K^2-lamdaBottom*(pmf_gauss(lamdaBottom)-pmf_gauss(lamdaTop));
noiseVar = noiseVar - pmf_gauss(lamdaTop)*(lamdaTop-lamdaBottom) + qfunc(lamdaTop)*(lamdaTop-lamdaBottom)^2;
noiseVar = signalVar*noiseVar;