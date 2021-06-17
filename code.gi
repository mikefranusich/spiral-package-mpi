ImportAll(simt);

Class(mpi_rcperm, call);

Class(MPICUDACodeGenMixin, rec(
    MPIRCPrm := (self,o,y,x,opts) >> simt_block(ApplyFunc(mpi_rcperm, [y,x]::o.func.params{[1..2]}::[Length(o.func.params[3])]::o.func.params[3]))
));


Class(MPICUDAUnparserMixin, rec(
    mpi_rcperm := (self, o, i, is) >> Print(Blanks(i), "fftx_mpi_rcperm", self.pinfix(o.args, ", "), ";\n")  
));


FixUpMPICUDAPerm := function(c)
    local mpicalls, mpinames, cucalls, mc, cc, ic;  

    mpicalls := Collect(c, @(1, specifiers_func, e->IsBound(e.cmd.cmds) and ObjId(e.cmd.cmds[1]) in [call, mpi_rcperm]));
    mpinames := List(mpicalls, i->i.id);
    cucalls := Collect(c, cu_call);
    
    for mc in mpicalls do
        SubstVars(mc, FoldR(Zip2(mc.params, Filtered(cucalls, cc->cc.func = mc.id)[1].args), (a,b) -> CopyFields(a, rec((b[1].id) := V(b[2]))), rec()));
    od;
    
    c:= SubstTopDown(c, @(1, cu_call, e->e.func in mpinames),
        e->Filtered(mpicalls, i->i.id = @(1).val.func)[1].cmd.cmds[1]
    );
    
    c := SubstTopDown(c, @(1, specifiers_func, e->IsBound(e.cmd.cmds) and ObjId(e.cmd.cmds[1]) in [call, mpi_rcperm]),
        e->skip());
        
    c := SubstBottomUp(c, @(1, chain), e-> chain(Filtered(@(1).val.cmds, c->ObjId(c) <> skip)));
    
    cc := Collect(c, @(1, call, e->IsBound(e.args[1].codegen)));
    ic := List(cc, _c -> _c.args[1].codegen.init());
    c := SubstBottomUp(c, @@(1, chain, (e,cx)->IsBound(cx.func) and Length(cx.func) = 1 and cx.func[1].id = "init"),
           e->chain(@@(1).val.cmds::ic));
        
    c := FoldL(List(cc, d->d.args[1].codegen), (a,b)->b.data(a), chain(c.cmds));
    c := program(decl(c.free(), c));
    c := SubstBottomUp(c, @(1, decl), 
        e -> decl(Filtered(@(1).val.vars, v -> v in @(1).val.cmd.free()), @(1).val.cmd));
        
    return c;    
end;
