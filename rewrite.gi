Class(RulesMPIPromoteNT, RuleSet);

RewriteRules(RulesMPIPromoteNT, rec(
    DAG_remove_mpi := Rule([@(1,TDAG), [@(2, TDAGNode), @(3), 
            @(4).cond(e->IsList(e) and Length(e) = 1 and ForAll(Collect(e[1], var), v-> v.id in ["Y", "Yptr"])), 
            @(5).cond(e->IsList(e) and Length(e) = 1 and ForAll(Collect(e[1], var), v-> v.id in ["X", "Xptr"])),
            ...],...],
        e-> @(2).val.params[1].withTags([AParMPI(rec(Y := @(4).val, X := @(5).val))])),
        
    TFunc_mergeTags := Rule([@(1, TFCall, e->e.hasTags()), @(2).cond(e->e.hasTags() and let(tg := e.getTags()[1], IsBound(tg.isParMPI) and tg.isParMPI)), ...], 
        e-> ApplyFunc(TFCall, [@(2).val.setTags([])]::Drop(@(1).val.params,1)).withTags(
            [ApplyFunc(ObjId(@(1).val.getTags()[1]), @(1).val.getTags()[1].params::@(2).val.getTags()[1].params)]::Drop(@(1).val.getTags(),1))
    )
));

RewriteRules(RulesRC, rec(
    RC_MPISum := Rule([RC, @(1, MPISum)],
        e -> let(s:=@(1).val, MPISum(s.nthreads, s.tid, s.var, s.domain, RC(s.child(1))))),
        
    RC_MPIPrm := Rule([RC, @(1,MPIPrm)], 
        e -> let(MPIRCPrm(@(1).val.func))),

    RC_MPIEmbed := Rule([RC, @(1,MPIPrm)], 
        e -> let(MPIRCPrm(@(1).val.func))),

    RC_MPIGath := Rule([RC, @(1, MPIGath)], e -> MPIGath(fTensor(@(1).val.func, fId(2)))),

    RC_MPIScat := Rule([RC, @(1, MPIScat)], e -> MPIScat(fTensor(@(1).val.func, fId(2))))
        
));



Class(LocalizeMPIRules, RuleSet);
RewriteRules(LocalizeMPIRules, rec(
    merge_GPS := ARule(Compose, [@(1, MPIGath), @(2, MPIRCPrm), @(3, MPIScat)],
        e->let(p := @(1).val.free()[1].range, f := @(2).val.func.params,
               [ MPIRCPrm(ApplyFunc(ObjId(@(2).val.func), [f[1]/p]::Drop(f, 1)))])),

    merge_GRCDS := ARule(Compose, [@(1, MPIGath), @(2, RCDiag), @(3, MPIScat)],
        e->let(p := @(1).val.free()[1].range, rcdf := @(2).val.element,
               [ RCDiag(FDataOfs(rcdf.var, rcdf.len/p, rcdf.ofs)) ])),


    drop_rcperm := ARule(Compose, [@(1, MPIGath), @(2, MPIRCPrm), @(3, MPIRCPrm), @(4, MPIScat)],
                         e-> let(
			          p := @(1).val.free()[1].range, f := @(2).val.func.params, l := f[3], param3 := [l[1]/p]::Drop(l, 1),
			         [ MPIRCPrm(ApplyFunc(ObjId(@(2).val.func), [f[1]/p]::DropLast(Drop(f, 1),1)::[param3]))])),



    drop_NoPull := Rule(@(1, NoPull), e->e.child(1)),
    
    drop_MPISum := Rule(@(1, MPISum), e->e.child(1))
    
));

Class(LocalizeMPIRulesCleanup, RuleSet);
RewriteRules(LocalizeMPIRulesCleanup, rec(
   drop_MPIGath_orScath := ARule(Compose, [@(1, [MPIGath, MPIScat])],
    			 e-> let( <#Error("ERROR HERE\n", e, "\n"),#> [])),
			 				  
));





FullLocalizeMPI := [MergedRuleSet(RulesSums, LocalizeMPIRules), LocalizeMPIRulesCleanup];



