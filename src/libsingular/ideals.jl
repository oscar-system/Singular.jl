
function id_Slimgb(I::ideal, R::ring; complete_reduction::Bool=false)
   if complete_reduction
      crbit = icxx"""Sy_bit(OPT_REDSB);"""
   else
      crbit = Cuint(0);
   end
   icxx"""ideal id = NULL;
          if (!idIs0($I))
          {
             intvec * n = NULL;
             tHomog h = testHomog;
             const ring origin = currRing;
             unsigned int save_opt = si_opt_1;
             si_opt_1 |= $crbit;
             rChangeCurrRing($R);
             id = t_rep_gb($R, $I, $I->rank);
             si_opt_1 = save_opt;
             rChangeCurrRing(origin);
             if (n != NULL)
                delete n;
          } else
             id = idInit(0, $I->rank);
          id;
       """
end

function id_Std(I::ideal, R::ring; complete_reduction::Bool=false)
   if complete_reduction
      crbit = icxx"""Sy_bit(OPT_REDSB);"""
   else
      crbit = Cuint(0);
   end
   icxx"""ideal id = NULL;
          if (!idIs0($I))
          {
             intvec * n = NULL;
             tHomog h = testHomog;
             const ring origin = currRing;
             unsigned int save_opt = si_opt_1;
             si_opt_1 |= $crbit;
             rChangeCurrRing($R);
             id = kStd($I, $R->qideal, h, &n);
             si_opt_1 = save_opt;
             rChangeCurrRing(origin);
             if (n != NULL)
                delete n;
          } else
             id = idInit(0, $I->rank);
          id;
       """
end



function id_fres(I::ideal, n::Cint, method::String, R::ring)
   s = icxx"""const ring origin = currRing;
         rChangeCurrRing($R);
         syStrategy s = syFrank($I, $n, $method);
         rChangeCurrRing(origin);
         s;
      """
   r = icxx"""$s->minres;"""
   minimal = true
   if r == C_NULL
      r = icxx"""$s->fullres;"""
      minimal = false
   end
   length = icxx"""$s->length;"""
   r, Int(length), minimal
end

function id_sres(I::ideal, n::Cint, R::ring)
   s = icxx"""const ring origin = currRing;
         rChangeCurrRing($R);
         syStrategy s = sySchreyer($I, $n);
         rChangeCurrRing(origin);
         s;
      """
   r = icxx"""$s->minres;"""
   minimal = true
   if r == C_NULL
      r = icxx"""$s->fullres;"""
      minimal = false
   end
   length = icxx"""$s->length;"""
   r, Int(length), minimal
end

function maGetPreimage(target::ring, map::ideal, id::ideal, source::ring)
   preimage_ptr = icxx"""
      sip_smap sing_map = { $map->m, (char *)"julia_ring", 1, $map->ncols };
      return maGetPreimage($target, &sing_map, $id, $source);
   """
   return preimage_ptr;
end
