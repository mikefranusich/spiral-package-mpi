_diceIt := (k,d) -> let(fct := Factors(k), l := Length(fct), lx:= Replicate(d-1, Int(l/d)), ll := l-Sum(lx), prt := lx::[ll],
    pfx := [0]::List([1..d], i -> Sum(prt{[1..i]})), prtl := List([1..d], i->Product(fct{[pfx[i]+1..pfx[i+1]]})), 
    Reversed(prtl));


NewRulesFor(MDDFT, rec(
    MDDFT_tSPL_pencil := rec(

        applicable := (self, t) >> Length(t.params[1]) > 1 and let(tg := t.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI),

        children := nt -> [let(tags := nt.getTags(), dd := Length(tags), peelof := Product(nt.getTags()[1].params[1]),
                              pencils := List([1..Length(nt.params)], d -> TMap(DFT(nt.params[1][d], nt.params[2]), 
                                  [peelof]::_diceIt(Product(nt.params[1]{[1..d-1]}::nt.params[1]{[d+1..Length(nt.params[1])]})/peelof, dd-1), APar, APar).withTags(nt.getTags())),
                              perms := List([1..Length(nt.params)], d -> TPrm(fCubeTranspose(Product(nt.params[1]), d, nt.params[1]))),
                              DropLast(Flat(Zip2(pencils, perms)), 1)
                           )],

        apply := (t, C, Nonterms) -> Compose(C),
        switch := true
    ),
    MDDFT_tSPL_pencil_cuFFT_MPI := rec(

        applicable := (self, t) >> Length(t.params[1]) > 1 and let(tg := t.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI),

        children := nt -> [[ TCompose([ 
                                TTensorII(DFT(nt.params[1][3], nt.params[2]), nt.params[1]{[1,2]}, ACubeRot_ZXY, ACubeRot_XYZ),
                                TPrm(fCubeTranspose(Product(nt.params[1]), "FFTX_MPI_3D_CUFFT_STAGE1", nt.params[1])), 
                                TTensorII(DFT(nt.params[1][2], nt.params[2]), nt.params[1]{[1,3]}, ACubeRot_ZXY, ACubeRot_XYZ),
                                TPrm(fCubeTranspose(Product(nt.params[1]), "FFTX_MPI_3D_CUFFT_STAGE2", nt.params[1])),
                                TTensorII(DFT(nt.params[1][1], nt.params[2]), nt.params[1]{[2,3]}, ACubeRot_ZXY, ACubeRot_XYZ)
                            ]).withTags(nt.getTags()) ]],

        apply := (t, C, Nonterms) -> C[1],
        switch := true
    )
));

NewRulesFor(TPrm, rec(
    TPrm_MPI := rec(
        forTransposition := false,
        children := nt -> [[  ]],
        apply := (nt, c, cnt) -> NoPull(MPIPrm(nt.params[1]))
    )
));

NewRulesFor(TTensorI, rec(
    IxA_MPI := rec(
        info := "IxA base",
        forTransposition := false,
        applicable := nt -> nt.hasTags() and let(tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI) and IsParPar(nt.params),
        children := nt -> [[ nt.params[1].withTags(Drop(nt.getTags(), 1)) ]],
        apply := (nt, c, cnt) ->  MPITensor(c[1], nt.params[2])
    )
));

NewRulesFor(TTensorII, rec(
    batch_cuFFT_cuberot_3D := rec(
        info := "IxA base",
        forTransposition := false,
        applicable := nt -> nt.hasTags() and let(tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI) and 
            nt.params[3] = ACubeRot_ZXY and nt.params[4] = ACubeRot_XYZ and
            ObjId(nt.params[1]) = DFT,
        children := nt -> [[ ]],
        apply := (nt, c, cnt) ->  let(K := nt.params[1].params[1], k := nt.params[1].params[2], M := nt.params[2][1], N := nt.params[2][2], 
                                      pg := nt.getTags()[1].params[1], p1 := pg[1], p2 := pg[2],
            MPITensor(
                CUFFTCall(
#                    Tensor(L(N*M/p1, M/p1), I(K/p2)) * Tensor(I(N), L(M*K/(p1*p2), M/p1)) * Tensor(I(p1*p2), DFT(K, k)), # this formula is only an approximation
                    TTensorII(nt.params[1], List([1..2], i-> nt.params[2][i]/nt.getTags()[1].params[1][i]), ACubeRot_ZXY, ACubeRot_XYZ),
                    rec(K := K, M := M, N := N, p1 := p1, p2 := p2, k := k)), 
                Product(pg)))
    )
));
