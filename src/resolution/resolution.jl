###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(r::sresolution) = r.base_ring

parent(r::sresolution) = SingularResolutionSet(r.base_ring)

function checkbounds(r::sresolution, i::Int)
   (i < 1 || i > r.len) && throw(BoundsError(I, i))
end

function getindex(r::sresolution, i::Int)
   checkbounds(r, i)
   R = base_ring(r)
   ptr = libSingular.getindex(r.ptr, Cint(i - 1))
   if ptr != C_NULL
      ptr = libSingular.id_Copy(ptr, R.ptr)
   end
   return Module(R, ptr)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::SingularResolutionSet)
   print(io, "Set of Singular Resolutions over ")
   show(io, R.base_ring)
end

function show(io::IO, r::sresolution)
   println(io, "Singular Resolution:")
   if r.len > 0
      print(io, "R^", rank(r[1]))
   end
   for i = 2:r.len
      m = r[i]
      if m.ptr == C_NULL
         break
      end
      print(io, " <- R^", rank(m))
   end
end
