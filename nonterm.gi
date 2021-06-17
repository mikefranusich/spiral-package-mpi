Class(ACubeRot_ZXY, ASingletonTag);
Class(ACubeRot_YZX, ASingletonTag);
Class(ACubeRot_XYZ, ASingletonTag);
Class(ACubeRot_ZYX, ASingletonTag);


Class(TTensorII, Tagged_tSPL_Container, rec(
    abbrevs :=  [ (nt, s, l, r) -> Checked(
        IsSPL(nt),
	[nt, s, l, r])
    ],

    dims := self >> self.params[1].dims()*Product(self.params[2]),

    SPLtSPL := (self, nt, P) >> Error("not implemented"),

    terminate := self >> Error("not implemented"),

    transpose := self >> let(p := self.params,
        TTensorII(p[1].transpose(), p[2], p[4], p[3]).withTags(self.getTags())),

    isReal := self >> self.params[1].isReal(),

    normalizedArithCost := self >>
        self.params[1].normalizedArithCost() * Product(self.params[2]),

    doNotMeasure := true,

    HashId := self >> let(
	p := self.params,
	h := When(IsBound(p[1].HashId), p[1].HashId(), p[1]),
        [h, p[2], p[3], p[4]] :: When(IsBound(self.tags), self.tags, [])
    ),
));

