Import(paradigms.smp);

Class(fCubeTranspose, PermClass, rec(
    domain := self >> self.params[1],
    range  := self >> self.params[1], 
    def    := (N, stage, cube) -> rec(),
    lambda := self >> let(i := Ind(self.params[1]), errExp(TInt)),
    transpose := self >> self,
    isIdentity := True
));


Class(fCubeEmbed, PermClass, rec(
    domain := self >> self.params[3][1],
    range  := self >> self.params[1],
    def    := (N, stage, cube) -> rec(),
    lambda := self >> let(i := Ind(self.params[1]), errExp(TInt)),
    transpose := self >> self,
    isIdentity := True
));

Class(MPIPrm, Prm,
 rec(
        dims := self >> let([self.func.range(), self.func.domain()])
));

Class(MPIRCPrm, Prm,
    rec(
       dims := self >> let(2*[self.func.range(), self.func.domain()])
));

Declare(MPITensor);

Class(MPISum, SMPSum);

Class(MPIGath, Gath);

Class(MPIScat, Scat);

Class(MPITensor, BaseMat, SumsBase, rec(
    dims := self >> self.L.dims() * self.P,
    isReal := self >> true,
    #-----------------------------------------------------------------------
    rChildren := self >> [self.L, self.P],
    rSetChild := rSetChildFields("L", "P"),
    #-----------------------------------------------------------------------
    new := (self, L, P) >> SPL(WithBases(self,
        rec(L   := L,
            P   := P,
            dimensions := L.dims()*P)
    )),

    #-----------------------------------------------------------------------
    transpose := self >> MPITensor(self.L.transpose(), self.P),
    #-----------------------------------------------------------------------

    print := (self,i,is) >> Print(self.name, "(",
        self.L.print(i+is,is), ", ", self.P, ")"),
    #-----------------------------------------------------------------------
    toAMat := self >> Tensor(I(self.P), self.L).toAMat(),
    #-----------------------------------------------------------------------
    isPermutation := False
));

Class(MPISumsgenMixin, rec(
    MPITensor := (self, o, opts) >> let(p:= o.P, i:= Ind(p),
        MPISum(p, i, p, 
        MPIScat(fTensor(fBase(i), fId(Rows(o.L)))) * o.L * MPIGath(fTensor(fBase(i), fId(Cols(o.L))))
    ))
));

