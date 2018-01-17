function sy_Delete(r::resolvente, i::Cint, R::ring)
   icxx"""id_Delete($r + $i, $R);"""
end

function getindex(r::resolvente, i::Cint)
   icxx"""(ideal) $r[$i];"""
end

