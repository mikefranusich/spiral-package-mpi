Load(fftx);
ImportAll(fftx);
Load(mpi);
Import(mpi);

# process grid and size
#pg := [4,4];
pg := [32,32];
d := 3;

#N := Replicate(d, 16);
N := Replicate(d, 1024);
procGrid := MPIProcGridND(pg);

name := "mdconv"::StringInt(d)::"d_"::
    StringInt(N[1])::ApplyFunc(ConcatenationString, List(Drop(N, 1), s->"x"::StringInt(s)))::"_mpi"::
    StringInt(pg[1])::ApplyFunc(ConcatenationString, List(Drop(pg, 1), s->"x"::StringInt(s)));
    
# machine configuration
conf := LocalConfig.mpi.confMPIGPUDevice(procGrid);

# on each node we have a local box of 1024 x 256 x 256 
Nlocal := N{[1]}::List([2, 3], i-> N[i]/pg[i-1]);
localBrick := TArrayND(TComplex, Nlocal, dimXYZ);

# thus our global box is 1024 x 1024 x 1024
dataLayout := TPtrGlobalArrayND(procGrid, localBrick);
Xglobal := tcast(dataLayout, X);
Yglobal := tcast(dataLayout, Y);

T1global := var("T1", dataLayout);
T2global := var("T2", dataLayout);
symvar := var("symvar", dataLayout);

# the MDDFT
mddft := TRC(MDDFT(N, -1));
imddft := TRC(MDDFT(N, 1));
#rcdiag := RCDiag(FDataOfs(symvar, 2*Product(N), 0).setRange(TReal));
rcdiag := RCDiag(FDataOfs(tcast(TPtr(TReal), symvar), 2*Product(N), 0));

# the transform
t := TFCall(
    TDecl(
        TDAG([
            TDAGNode(mddft, T1global, Xglobal),
            TDAGNode(rcdiag, T2global, T1global),
            TDAGNode(imddft, Yglobal, T2global)
        ]), [ T1global, T2global ]), 
        rec(fname := name, params := [ symvar ])
);

# derive the opts and tag the transform
opts := conf.getOpts(t);
tt := opts.tagIt(t);
c := opts.fftxGen(tt);

opts.prettyPrint(c);
PrintTo(name::".cu", opts.prettyPrint(c));
                                            
