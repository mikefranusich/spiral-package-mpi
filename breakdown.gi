_diceIt := (k,d) -> let(fct := Factors(k), l := Length(fct), lx:= Replicate(d-1, Int(l/d)), ll := l-Sum(lx), prt := lx::[ll],
    pfx := [0]::List([1..d], i -> Sum(prt{[1..i]})), prtl := List([1..d], i->Product(fct{[pfx[i]+1..pfx[i+1]]})), 
    Reversed(prtl));


NewRulesFor(MDDFT, rec(

    MDDFT_tSPL_pencil := rec(

        applicable := (self, t) >> Length(t.params[1]) > 1 and let(tg := t.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI),

        children := nt -> [let(tags := nt.getTags(), dd := Length(tags), peelof := Product(nt.getTags()[1].params[1]),
		       	      Debug(),
		       	      #create pencils compute
		       	      pencils := List([1..Length(nt.params)], d -> TMap(DFT(nt.params[1][d], nt.params[2]), 
                                  [peelof]::_diceIt(Product(nt.params[1]{[1..d-1]}::nt.params[1]{[d+1..Length(nt.params[1])]})/peelof, dd-1), APar, APar).withTags(nt.getTags())),

			      #create embeddings
#			      embs := List([1..Length(nt.params)], d -> TPrm(fCubeEmbed(Product(nt.params[1]), "\"FFTX_MPI_EMBED_"::String(d)::"\"", [Product(),true])));

			      #create permutations
                              perms := List([1..Length(nt.params)], d -> TPrm(fCubeTranspose(Product(nt.params[1]), "\"FFTX_MPI_3D_CUFFT_STAGE"::String(d)::"\"", nt.params[1]))),


			      Print(pencils[1]), Error("--"),
                              DropLast(Flat(Zip2(pencils, perms)), 1)
                           )],

        apply := (t, C, Nonterms) -> Compose(C),
        switch := true
    ),
    
    MDDFT_tSPL_pencil_cuFFT_MPI := rec(

        applicable := (self, t) >> Length(t.params[1]) > 1 and let(tg := t.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI),

        children := nt -> let(Print("DFT Size: ",nt.params[1][3], " Batch:", nt.params[2], " ",nt.params[1]," ",nt.getTags(), "\n"),
		       	  [[ TCompose([ 
                                TTensorII(DFT(nt.params[1][3], nt.params[2]), nt.params[1]{[1,2]}, ACubeRot_ZXY, ACubeRot_XYZ),
                                TPrm(fCubeTranspose(Product(nt.params[1]), "FFTX_MPI_3D_CUFFT_STAGE1", nt.params[1])), 
                                TTensorII(DFT(nt.params[1][2], nt.params[2]), nt.params[1]{[3,1]}, ACubeRot_ZXY, ACubeRot_XYZ),
                                TPrm(fCubeTranspose(Product(nt.params[1]), "FFTX_MPI_3D_CUFFT_STAGE2", nt.params[1])),
                                TTensorII(DFT(nt.params[1][1], nt.params[2]), nt.params[1]{[2,3]}, ACubeRot_ZXY, ACubeRot_XYZ)
                            ]).withTags(nt.getTags()) ]]),

        apply := (t, C, Nonterms) -> C[1],
        switch := true
    )
));

NewRulesFor(TPrm, rec(
    TPrm_MPI := rec(
        forTransposition := false,
        children := nt -> [[  ]],
        apply := (nt, c, cnt) -> let(NoPull(MPIPrm(nt.params[1])))
    )
));

NewRulesFor(TTensorI, rec(
    IxA_MPI := rec(
        info := "IxA base",
        forTransposition := false,
        applicable := nt -> nt.hasTags() and let(tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI) and IsParPar(nt.params),
        children := nt -> [[ nt.params[1].withTags(Drop(nt.getTags(), 1)) ]],
        apply := (nt, c, cnt) ->  let(MPITensor(c[1], nt.params[2]))
    ),
    AxI_MPI := rec(
        info := "AxI base",	
     	forTransposition := false,
	applicable := nt -> nt.hasTags() and let (tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI) and IsVecVec(nt.params),
	children := nt -> [[ nt.params[1].withTags(Drop(nt.getTags(), 1)) ]],
	apply := (nt, c, cnt) -> MPITensor(c[1], nt.params[2])
     ),
 
    L_IxA_MPI := rec(
        forTransposition := false,
        
        # these config parameters need to be moved into the opts...
        mem := 1024*96,
        mem_per_pt := 2*8*2,
        max_threads := 2048,
        max_kernel := 18 * 18,
        _peelof := (self,n,m) >> Maximum(Filtered(self.mem_per_pt * Filtered(n*DivisorsInt(m), e-> e<self.max_threads), 
            f -> f < When(n >= self.max_kernel, self.mem/2, self.mem)))/(self.mem_per_pt*n),
        
        applicable := nt -> nt.hasTags() and let (tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI)  and IsVecPar(nt.params), # and nt.params[2] > 1,
        children := (self, nt) >> let(Error("L_IXA"), n := Rows(nt.params[1]), m:= nt.params[2], peelof := self._peelof(n,m), remainder := m/peelof,
            [[ TCompose([TL(Rows(nt)/peelof, n, 1, peelof), 
                TTensorI(
                    TCompose([TL(Rows(nt.params[1]) * peelof, Rows(nt.params[1])), TTensorI(nt.params[1], peelof, APar, APar)]),
                    remainder, APar, APar)]).withTags(nt.getTags()) ]]),
        apply := (nt, c, cnt) -> MPITensor(c[1], nt.params[2])
    ),

    LL_IxA_MPI := rec(
        forTransposition := false,
        
        # these config parameters need to be moved into the opts...
        mem := 1024*96,
        mem_per_pt := 2*8*2,
        max_threads := 2048,
        max_kernel := 18 * 18,
        _peelof := (self,n,m) >> Maximum(Filtered(self.mem_per_pt * Filtered(n*DivisorsInt(m), e-> e<self.max_threads), 
            f -> f < When(n >= self.max_kernel, self.mem/2, self.mem)))/(self.mem_per_pt*n),
        
        applicable := nt -> nt.hasTags() and let (tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI)  and IsParVec(nt.params), # and nt.params[2] > 1,
        children := (self, nt) >> let(
	    n := Rows(nt.params[1]), m:= nt.params[2], peelof := self._peelof(n,m), remainder := m/peelof,
            [[ TCompose([TL(Rows(nt)/peelof, n, 1, peelof), 
                TTensorI(
                    TCompose([TL(Rows(nt.params[1]) * peelof, Rows(nt.params[1])), TTensorI(nt.params[1], peelof, APar, APar)]),
                    remainder, APar, APar)]).withTags(nt.getTags()) ]]),
        apply := (nt, c, cnt) -> MPITensor(c[1], nt.params[2])
    )
));

NewRulesFor(TTensorII, rec(

batch_cuFFT_cuberot_3D := rec(
        info := "IxA base",
        forTransposition := false,
        applicable := nt -> nt.hasTags() and let(tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI)
            and nt.params[3]<>nt.params[4] #ACubeRot_ZXY and nt.params[4]=ACubeRot_XYZ
	    and  ObjId(nt.params[1]) = DFT,
        children := nt -> [[ ]],
        apply := (nt, c, cnt) ->  let(
				     #K := nt.params[1].params[1], k := nt.params[1].params[2], M := nt.params[2][1], N := nt.params[2][2],
				     K := nt.params[2][2], k := nt.params[1].params[2], M := nt.params[2][1], N := nt.params[1].params[1],
                                     pg := nt.getTags()[1].params[1], p1 := pg[1], p2 := pg[2],
				     indist:= K,
				     instride := 1, #M*K/Product(pg), 
				     outstride:= M*K/Product(pg), 
				     outdist := 1,
            MPITensor(
                CUFFTCall(
                    TTensorII(nt.params[1], List([1..2], i-> nt.params[2][i]/nt.getTags()[1].params[1][i]), nt.params[3], nt.params[4]),
                    rec(K := K, M := M, N := N, p1 := p2, p2 := p2, k := k, instride := instride, outstride := outstride, indist := indist, outdist := outdist)), 
                Product(pg)))
   ),

batch_cuFFT_cubenorot_3D := rec(
        info := "IxA base",
        forTransposition := false,
        applicable := nt -> nt.hasTags() and let(tg := nt.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI)
	    and nt.params[3]=nt.params[4]
	    and  ObjId(nt.params[1]) = DFT,
        children := nt -> [[ ]],
        apply := (nt, c, cnt) ->  let(
				     K := nt.params[2][2], k := nt.params[1].params[2], M := nt.params[2][1], N := nt.params[1].params[1],
                                     pg := nt.getTags()[1].params[1], p1 := pg[1], p2 := pg[2],
				     indist:= K,
				     instride := 1, #M*K/Product(pg), ##Product(nt.params[2]/Product(pg)), #N
				     outstride:= 1, #M*K/Product(pg), #nt.params[1].params[1],
				     outdist := K,
            MPITensor(
                CUFFTCall(
                    TTensorII(nt.params[1], List([1..2], i-> nt.params[2][i]/nt.getTags()[1].params[1][i]), nt.params[3], nt.params[4]),
                    rec(K := K, M := M, N := N, p1 := p2, p2 := p2, k := k, instride := instride, outstride := outstride, indist := indist, outdist := outdist)), 
                Product(pg)))
   ),
    

));
