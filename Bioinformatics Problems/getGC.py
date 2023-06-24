def getGC(strlist):
    counts = [(curStr.count("C")+curStr.count("G"),len(curStr)) for curStr in strlist]
    res = []
    for (gc,tot) in counts:
        gc = float(gc)
        tot = float(tot)
        if tot != 0:
            res.append(gc/tot)
    if (res):
        maxCount = max(res)
        indexofMax = res.index(maxCount)
        strMax = strlist[indexofMax]
        return (maxCount,strMax,indexofMax)
    return (0.0,"",0)

def main():
    s = "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"
    s1 = "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"
    s2 = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
    li = [s,s1,s2]
    res = getGC(li)
    print res

