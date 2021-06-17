Load(fftx);
ImportAll(fftx);
Load(mpi);
Import(mpi);

# process grid and size
pg := [2,2];
#pg := [32,32];
d := 3;
useCUFFT := true;

N := Replicate(d, 32);
#N := Replicate(d, 1024);
procGrid := MPIProcGridND(pg);
copts := rec(useCUFFT := useCUFFT);

name := "mddft"::StringInt(d)::"d_"::
    StringInt(N[1])::ApplyFunc(ConcatenationString, List(Drop(N, 1), s->"x"::StringInt(s)))::"_mpi"::
    StringInt(pg[1])::ApplyFunc(ConcatenationString, List(Drop(pg, 1), s->"x"::StringInt(s)))::When(useCUFFT, "_cufft", "");

# machine configuration
#conf := LocalConfig.mpi.confMPIGPUDevice(procGrid);
conf := LocalConfig.mpi.confMPIGPUDevice(procGrid, copts);

# on each node we have a local box of 1024 x 256 x 256 
Nlocal := N{[1]}::List([2, 3], i-> N[i]/pg[i-1]);
localBrick := TArrayND(TComplex, Nlocal, dimXYZ);

# thus our global box is 1024 x 1024 x 1024
dataLayout := TPtrGlobalArrayND(procGrid, localBrick);
Xglobal := tcast(dataLayout, X);
Yglobal := tcast(dataLayout, Y);

# the MDDFT
mddft := TRC(MDDFT(N, -1));

# the transform
t := TFCall(
    TDAG([TDAGNode(mddft, Yglobal, Xglobal)]),
    rec(fname := name, params := [ ])
);

# derive the opts and tag the transform
opts := conf.getOpts(t);
tt := opts.tagIt(t);
c := opts.fftxGen(tt);

# Temporary infrastructure for end-to-end demo 
init_comm_fn_name := var("init_2d_comms"); 
lst := [init_comm_fn_name]::pg::N; 
init_comm := ApplyFunc(call, lst); 
cc := Collect(c, func); 
cc[1].cmd := chain(init_comm, cc[1].cmd);
destroy_comm_fn_name := var("destroy_2d_comms"); 
destroy_comm := ApplyFunc(call, [destroy_comm_fn_name]); 
cc[3].cmd := chain(destroy_comm, cc[3].cmd);

cx := program(cc);
cx.ruletree := c.ruletree;
opts.prettyPrint(cx);
PrintTo(name::".cu", opts.prettyPrint(cx));
                                            
