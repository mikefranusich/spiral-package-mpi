MPICommunicator := grd -> AParMPI(grd);

Class(TArrayND, TArrayBase, rec(
    isArrayT := true,
    toPtrType := self >> TPtr(self.t),
    doHashValues := true,

     __call__ := (self, t, sizes, linearization) >>
        WithBases(self, rec(
        t    := Checked(IsType(t), t),
        sizes := sizes,
        linearization := linearization,
        operations := TypOps)),
    print := self >> Print(self.__name__, "(", self.t, ", ", self.sizes, ", ", self.linearization, ")"),
    
    rChildren := self >> [self.t, self.sizes, self.linearization],
    rSetChild := rSetChildFields("t", "sizes", "linearization"),
    free := self >> []
));

Class(TGlobalArrayND, TArrayBase, rec(
    isArrayT := true,
    toPtrType := self >> TPtr(self.t),
    doHashValues := true,

     __call__ := (self, pgrid, localArray) >>
        WithBases(self, rec(
        t    := localArray.t,
        pgrid := pgrid,
        localArray := localArray,
        operations := TypOps)),
    print := self >> Print(self.__name__, "(", self.pgrid, ", ", self.localArray, ")"),

    rChildren := self >> [self.t, self.pgrid, self.localArray],
    rSetChild := rSetChildFields("t", "pgrid", "localArray"),
    free := self >> []
));

Class(TPtrGlobalArrayND, TPtr, rec(
    doHashValues := true,

     __call__ := (self, pgrid, localArray) >>
        WithBases(self, rec(
        t    := localArray.t,
        pgrid := pgrid,
        localArray := localArray,
        operations := TypOps)),
    print := self >> Print(self.__name__, "(", self.pgrid, ", ", self.localArray, ")"),

    rChildren := self >> [self.t, self.pgrid, self.localArray],
    rSetChild := rSetChildFields("t", "pgrid", "localArray"),
    free := self >> []
));
