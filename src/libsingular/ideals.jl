function id_Slimgb(I::ideal, R::ring; complete_reduction::Bool=false)
    id_Slimgb_helper(I, R,complete_reduction)
end

function maGetPreimage(target::ring, map::ideal, id::ideal, source::ring)
   preimage_ptr = icxx"""
      sip_smap sing_map = { $map->m, (char *)"julia_ring", 1, $map->ncols };
      return maGetPreimage($target, &sing_map, $id, $source);
   """
   return preimage_ptr;
end
