# MPI Baseline options

ImportAll(fftx);
ImportAll(simt);

Declare(MPIOpts);
Declare(MPIGPUDeviceOpts);

Class(MPIDefaultConf, LocalConfig, rec(
    getOpts := meth(self, t)
                   local _opts;
                   _opts := CopyFields(MPIOpts, rec(tags := self.tags));
                   _opts.breakdownRules.TTensorI := [IxA_base];
                   _opts.breakdownRules.TRC := [CopyFields(TRC_tag, rec(applicable := True))];
                   _opts.breakdownRules.TPrm := [TPrm_MPI];
                   return _opts;    
               end,
                   
    operations := rec(Print := s -> Print("<MPI Default Configuration>")),
    tags := []
));

    new_decl := meth(self,o,i,is)
            local arrays, other, l, arri, myMem;
   	     [arrays, other] := SplitBy(o.vars, x->IsArray(x.t));
	     DoForAll(arrays, v -> Print(Blanks(i),
             When(self.opts.arrayBufModifier <> "", self.opts.arrayBufModifier::" ", ""),
                                 self.declare(v.t, v, i, is), ";\n"));

             if (Length(other)>0) then
	         other:=SortRecordList(other,x->x.t);
	         for l in other do
                     Sort(l, (a,b)->a.id < b.id);
                     Print(Blanks(i),
                           When(self.opts.arrayBufModifier <> "", self.opts.arrayBufModifier::" ", ""),
		   	   self.declare(l[1].t, l, i, is), ";\n");										          od;																		fi;

             self(o.cmd, i, is);

	     end;




Class(MPIGPUDeviceConf, fftx.platforms.cuda.FFTXCUDADeviceOpts, rec(
    getOpts := meth(self, t)
                   local opts;
                   
                   opts := batchFftCUDADeviceOpts();
                   opts.operations := rec(Print := s -> Print("<MPI Distributed FFT/CUDA Device options>"));

                   opts.tags := self.tags::opts.tags;
                   opts.breakdownRules.TTensorI := [IxA_MPI, IxA_base, AxI_MPI, L_IxA_MPI, LL_IxA_MPI, fftx.platforms.cuda.L_IxA_SIMT, fftx.platforms.cuda.IxA_SIMT]::opts.breakdownRules.TTensorI;
                   opts.breakdownRules.TTensorII := [batch_cuFFT_cuberot_3D, batch_cuFFT_cubenorot_3D];
                   opts.breakdownRules.TRC := [CopyFields(TRC_tag, rec(applicable := True))];
                   opts.breakdownRules.TPrm := [TPrm_MPI];
                   
#                   opts.breakdownRules.MDDFT := When(Length(self.copts) > 0 and IsRec(self.copts[1]) and IsBound(self.copts[1].useCUFFT) and self.copts[1].useCUFFT, 
#                       [ MDDFT_tSPL_pencil_cuFFT_MPI ], [ MDDFT_tSPL_pencil ]);

	           opts.includes := [ "\"fftx_mpi.hpp\"", "<stdint.h>" ],
                   opts.preProcess := (self, t) >> RulesMPIPromoteNT(t);
                   
                   opts.preProcess := t -> ApplyStrategy(t, 
                    [ RulesMPIPromoteNT,
                      RulesFFTXPromoteNT, 
                      RulesFFTXPromoteNT_Cleanup ],
                       BUA, opts);                   
                   
                   opts.fftxGen := MPIGPUDeviceOpts.fftxGen;
                   opts._sumsRuleTree := opts.sumsRuleTree;
                   opts.sumsRuleTree := MPIGPUDeviceOpts.sumsRuleTree;
                   
                   opts.codegen.MPIRCPrm := MPICUDACodeGenMixin.MPIRCPrm;
                   opts.unparser.mpi_rcperm := MPICUDAUnparserMixin.mpi_rcperm;
                   opts.sumsgen.MPITensor := MPISumsgenMixin.MPITensor;
                   
                   opts.dynamicDeviceMemory := true;

#		   opts.codegen.arrayDataModifier := "static";
#		   opts.codegen.arrayBufModifier := "static";

		   opts.arrayBufModifier := "static";
		   opts.arrayDataModifier := "static";



#		   arrayBufModifier := "static";
#		   arrayDataModifier := "static";

		   opts.unparser.data := (self,o,i,is) >> Print(
		                                   When(not IsArrayT(o.var.t), Print("static", Blanks(i), self.genData(o.var, o.value))),
 			                            self(o.cmd, i, is));
		   opts.unparser.decl := new_decl;				 
<#				       
#>

                   return opts;    
               end,
                   
    operations := rec(Print := s -> Print("<MPI GPU Device Configuration>")),
    tags := []
));


Class(MPIGPUDeviceOpts, SpiralDefaults, rec(
    tags := [],
    includes := [ "\"fftx_mpi.hpp\"", "<stdint.h>" ],
    tagIt := (self, t) >> t.withTags(self.tags),
    operations := rec(Print := s -> Print("<MPI default options>")),
        
    preProcess := (self, t) >> RulesMPIPromoteNT(t),
    
    search := (self, t) >> RandomRuleTree(t, self),
    
    sumsRuleTree := meth(self, rt) 
        local ss;
        ss := self._sumsRuleTree(rt);
        ss := ApplyStrategy(ss, FullLocalizeMPI, BUA, self);
        return ss;
    end,    
    
    fftxGen := meth(self, tt)
        local _tt, rt, ss, c, plist, plist1, tf, tf2;
    
        _tt := self.preProcess(tt);
        rt := self.search(_tt);
        ss := self.sumsRuleTree(rt);

        self.symbol := When(IsBound(rt.node.params[2].params), rt.node.params[2].params, []);
        
        c := CodeSums(ss, self);
        c := FixUpMPICUDAPerm(c);
        c := fixReplicatedData(c, self);
        
        plist := Collect(c, @(1,func, e->e.id="transform"))[1].params;
        plist1 := _orderedUniquify(plist);
        if Length(plist1) = 2 then 
            plist1 := plist1 :: self.symbol; 
        fi;
        c := SubstTopDown(c, @(1).cond(e->IsList(e) and e = plist), e->plist1);
        
        if self.symbol <> [] then
            c := SubstTopDown(c, @(1,func, e->e.id="transform"),
                e->SubstVars(e, rec((self.symbol[1].id) := var("_"::self.symbol[1].id, TPtr(TReal)))));
    
            c := SubstTopDown(c, @(1,specifiers_func),
                e->SubstVars(e, rec((self.symbol[1].id) := var("__"::self.symbol[1].id, TPtr(TReal)))));
    
            c := SubstTopDown(c, @(1, decl, e-> self.symbol[1] in e.vars),
                e-> decl(Filtered(e.vars, j->j <> self.symbol[1]), e.cmd));
        fi;
        
        c.ruletree := rt;
        c.ss := ss;
        return c;
    end
));


buildMPIGPUDeviceConf := function(arg)
    local conf;

    conf := CopyFields(MPIGPUDeviceConf, 
        rec(
            tags := arg{[1]}, 
            copts := Drop(arg,1),
            gpuconf := fftx.platforms.cuda.FFTXCUDADeviceDefaultConf
        ));
    
    return conf;
end;


Class(MPIOpts, SpiralDefaults, rec(
    tags := [],
    includes := [ "\"fftx_mpi.hpp\"", "<stdint.h>" ],
    tagIt := (self, t) >> t.withTags(self.tags),
        operations := rec(Print := s -> Print("<MPI default options>")),
        
    preProcess := (self, t) >> RulesMPIPromoteNT(t),
    
    search := (self, t) >> RandomRuleTree(t, self),
));


Class(MPIGlobals, rec(
    confMPI := (arg) >> CopyFields(MPIDefaultConf, rec(tags := [arg[2]])),
    defaultConf := arg >> ApplyFunc(arg[1].confMPI, Drop(arg, 1)),
    confMPIGPUDevice := arg >> ApplyFunc(buildMPIGPUDeviceConf, Drop(arg, 1))
));

spiral.LocalConfig.mpi := MPIGlobals;


