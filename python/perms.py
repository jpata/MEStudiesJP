import ROOT, sys, numpy

tf = ROOT.TFile(sys.argv[1])

tt = tf.Get("tree")

ps = numpy.zeros(60, dtype=numpy.int32)
n = numpy.zeros(1, dtype=numpy.int32)
tt.SetBranchAddress("perm_to_gen_s", ps)
tt.SetBranchAddress("nPermut_s", n)


for i in range(tt.GetEntries()):
    
    tt.GetEntry(i)
    
    s = ""

    for j in range(n[0]):
        s += str(ps[j]) + ","
    if len(s)>0:
        print s