## NOT USED

import numpy as np
import pandas as pd
from statsmodels.robust import mad
import scipy

def rep(i, x):
    return [i]*x

def pcf(df):
    ## takes in df[df["chromosome"] == c][["chromosome", "pos", s]]
    
    def filterMarkS4(x, kmin, L, L2, frac1, frac2, frac3, thres):
        lengdeArr = len(x)
        xc = np.cumsum(x)
        xc = [0, xc]
        ind11 = list(range(1, lengdeArr - 6 * L + 2))
        ind12 = [i+L for i in ind11] #ind11 + L
        ind13 = [i+(3*L) for i in ind11] #ind11 + 3 * L
        ind14 = [i+(5*L) for i in ind11]#ind11 + 5 * L
        ind15 = [i+(6*L) for i in ind11]#ind11 + 6 * L
        cost1 = [abs(4 * c-a-b-d-e) for a,b,c,d,e in zip(ind11, ind12, ind13,ind14,ind15)]
        cost1 = pd.Series([rep(0, 3 * L - 1)] + cost1 + [rep(0, 3 * L)])
        print(cost1)
        ci = len(cost1)
        in1 = list(range(1, (lengdeArr - 6)+1))
        in2 = [i + 1 for i in in1]
        in3 = [i + 2 for i in in1]
        in4 = [i + 3 for i in in1]
        in5 = [i + 4 for i in in1]
        in6 = [i + 5 for i in in1]
        in7 = [i + 6 for i in in1]
        test = np.maximum(cost1, cost1[in7])
        test = list(rep(0, 3), test, rep(0, 3))
        cost1B = cost1[cost1 >= thres * test]
        frac1B = min(0.8, frac1 * len(cost1)/len(cost1B))
        limit = quantile(cost1B, (1 - frac1B), names = False)
        mark = (cost1 > limit) & (cost1 > 0.9 * test)
        ind21 = list(range(1, lengdeArr - 6 * L2 + 2))
        ind22 = ind21 + L2
        ind23 = ind21 + 3 * L2
        ind24 = ind21 + 5 * L2
        ind25 = ind21 + 6 * L2
        cost2 = abs(4 * xc[ind23] - xc[ind21] - xc[ind22] - xc[ind24] - 
            xc[ind25])
        limit2 = quantile(cost2, (1 - frac2), names = False)
        mark2 = (cost2 > limit2)
        mark2 = list(rep(0, 3 * L2 - 1), mark2, rep(0, 3 * L2))
        if (3 * L > kmin):
            mark[kmin:(3 * L - 1)] = True
            mark[(lengdeArr - 3 * L + 1):(lengdeArr - kmin)] = True
        else:
            mark[kmin] = True
            mark[lengdeArr - kmin] = True
        if kmin > 1:
            ind1 = list(range(1, (lengdeArr - 3 * kmin + 1)+1))
            ind2 = ind1 + 3 * kmin
            ind3 = ind1 + kmin
            ind4 = ind1 + 2 * kmin
            shortAb = abs(3 * (xc[ind4] - xc[ind3]) - (xc[ind2] - 
                xc[ind1]))
            in1 = list(range(1, (len(shortAb) - 6)+1))
            in2 = in1 + 1
            in3 = in1 + 2
            in4 = in1 + 3
            in5 = in1 + 4
            in6 = in1 + 5
            in7 = in1 + 6
            test = pmax(shortAb[in1], shortAb[in2], shortAb[in3], 
                shortAb[in4], shortAb[in5], shortAb[in6], shortAb[in7])
            test = list(rep(0, 3), test, rep(0, 3))
            cost1C = shortAb[shortAb >= thres * test]
            frac1C = min(0.8, frac3 * len(shortAb)/len(cost1C))
            limit3 = quantile(cost1C, (1 - frac1C), names = False)
            markH1 = (shortAb > limit3) & (shortAb > thres * test)
            markH2 = list(rep(False, (kmin - 1)), markH1, rep(False, 
                2 * kmin))
            markH3 = list(rep(False, (2 * kmin - 1)), markH1, rep(False, 
                kmin))
            mark = mark | mark2 | markH2 | markH3
        else:
            mark = mark | mark2
        if (3 * L) > kmin:
            mark[1:(kmin - 1)] = False
            mark[kmin:(3 * L - 1)] = True
            mark[(lengdeArr - 3 * L + 1):(lengdeArr - kmin)] = True
            mark[(lengdeArr - kmin + 1):(lengdeArr - 1)] = False
            mark[lengdeArr] = True
        else:
            mark[1:(kmin - 1)] = False
            mark[(lengdeArr - kmin + 1):(lengdeArr - 1)] = False
            mark[lengdeArr] = True
            mark[kmin] = True
            mark[lengdeArr - kmin] = True
        
        return mark

    def findEst(bestSplit, N, Nr, Sum, yest):
        n = N
        lengde = rep(0, N)
        antInt = 0
        while (n > 0):
            antInt = antInt + 1
            lengde[antInt] = n - bestSplit[n]
            n = bestSplit[n]
        lengde = lengde[antInt:1]
        lengdeOrig = rep(0, antInt)
        startOrig = rep(1, antInt + 1)
        verdi = rep(0, antInt)
        start = rep(1, antInt + 1)
        for i in range(1, antInt+1):
            start[i + 1] = start[i] + lengde[i]
            lengdeOrig[i] = sum(Nr[start[i]:(start[i + 1] - 1)])
            startOrig[i + 1] = startOrig[i] + lengdeOrig[i]
            verdi[i] = sum(Sum[start[i]:(start[i + 1] - 1)])/lengdeOrig[i]
        if yest:
            yhat = rep(0, startOrig[antInt + 1] - 1)
            for i in range(1, antInt+1):
                yhat[startOrig[i]:(startOrig[i + 1] - 1)] = verdi[i]
            startOrig = startOrig[1:antInt]
            return {"Lengde": lengdeOrig, "sta": startOrig, "mean": verdi,
                    "nIntervals": antInt, "yhat": yhat}
        else:
            startOrig = startOrig[1:antInt]
            return {"Lengde": lengdeOrig, "sta": startOrig, "mean": verdi,
                    "nIntervals": antInt}
        

    def PottsCompact(kmin, gamma, nr, res, yest):
            N = len(nr)
            Ant = rep(0, N)
            Sum = rep(0, N)
            Cost = rep(0, N)
            if (sum(nr) < 2 * kmin):
                estim = sum(res)/sum(nr)
                return estim
            initAnt = nr[1]
            initSum = res[1]
            initAve = initSum/initAnt
            bestCost = rep(0, N)
            bestCost[1] = (-initSum * initAve)
            bestSplit = rep(0, N)
            k = 2
            while (sum(nr[1:k]) < 2 * kmin):
                Ant[2:k] = Ant[2:k] + nr[k]
                Sum[2:k] = Sum[2:k] + res[k]
                bestCost[k] = (-(initSum + Sum[2])^2/(initAnt + Ant[2]))
                k <- k + 1
            for n in range(k, N+1):
                Ant[2:n] = Ant[2:n] + nr[n]
                Sum[2:n] = Sum[2:n] + res[n]
                limit = n
                while (limit > 2 & Ant[limit] < kmin):
                    limit = limit - 1
                Cost[2:limit] = bestCost[1:limit - 1] - Sum[2:limit]^2/Ant[2:limit]
                Pos = which.min(Cost[2:limit]) + 1
                cost = Cost[Pos] + gamma
                totCost = (-(Sum[2] + initSum)^2/(Ant[2] + initAnt))
                if (totCost < cost):
                    Pos = 1
                    cost = totCost
                bestCost[n] = cost
                bestSplit[n] = Pos - 1
            
            if yest:
                yhat = rep(0, N)
                res = findEst(bestSplit, N, nr, res, True)
            else:
                res = findEst(bestSplit, N, nr, res, False)
            return res

    def compact(y, mark):
        tell = list(range(1, len(y)+1))
        cCTell = tell[mark]
        Ncomp = len(cCTell)
        lowTell = list(0, cCTell[1:(Ncomp - 1)])
        ant = cCTell - lowTell
        cy = cumsum(y)
        cCcy = cy[mark]
        lowcy = list(0, cCcy[1:(Ncomp - 1)])
        sumcy = cCcy - lowcy
        return {"Nr": ant, "Sum": sumcy}

    def runFastPcf(x, kmin, gamma, fract1, fract2, yest=True):
        antGen = len(x)
        mark = filterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9)
        mark[antGen] = TRUE
        dense = compact(y, mark)
        result = PottsCompact(kmin, gamma, dense["Nr"], dense["Sum"], yest)
        return result 

    def runPcfSubset(x, kmin, gamma, frac1, frac2, yest):
        SUBSIZE = 5000
        antGen = len(x)
        mark = filterMarkS4(x, kmin, 8, 1, frac1, frac2, 0.02, 0.9)
        markInit = list(mark[1:(SUBSIZE - 1)], True)
        compX = compact(x[1:SUBSIZE], markInit)
        mark2 = rep(False, antGen)
        mark2[1:SUBSIZE] = markWithPotts(kmin, gamma, compX["Nr"], 
            compX["Sum"], SUBSIZE)
        mark2[4 * SUBSIZE/5] = True
        start = 4 * SUBSIZE/5 + 1
        while (start + SUBSIZE < antGen):
            slutt = start + SUBSIZE - 1
            markSub = list(mark2[1:(start - 1)], mark[start:slutt])
            markSub[slutt] = True
            compX = compact(x[1:slutt], markSub)
            mark2[1:slutt] = markWithPotts(kmin, gamma, compX["Nr"], 
                compX["Sum"], slutt)
            start = start + 4 * SUBSIZE/5
            mark2[start - 1] = True
        
        markSub = list(mark2[1:(start - 1)], mark[start:antGen])
        compX = compact(x, markSub)
        result = PottsCompact(kmin, gamma, compX["Nr"], compX["Sum"], 
            yest)
        return result
    
    cols = list(df.columns)
    sample = cols[-1]
    nProbe = len(df["pos"])
    gamma = 4000
    segments = pd.DataFrame(columns=["sampleID", "chrom", "arm", "start", "end", "nprobes", "mean"])
    sd = scipy.stats.median_absolute_deviation(df[sample])
    print(sd)
    x = [abs(i-sd) for i in df[sample]]

    xLength = len(x)
    kmin = 5
    yest = True
    if (xLength < 1000):
        result = runFastPcf(x, kmin, gamma, 0.15, 0.15, yest)
    else:
        if (xLength < 15000):
            result = runFastPcf(x, kmin, gamma, 0.12, 0.05, yest)
        else:
            result = runPcfSubset(x, kmin, gamma, 0.12, 0.05, yest)
    return result

