function sy_Delete(r::resolvente, i::Cint, R::ring)
   icxx"""id_Delete($r + $i, $R);"""
end

function getindex(r::resolvente, i::Cint)
   icxx"""(ideal) $r[$i];"""
end

function syMinimize(r::resolvente, length::Cint, R::ring)
   icxx"""const ring origin = currRing;
          syStrategy temp = (syStrategy) omAlloc0(sizeof(ssyStrategy));
          resolvente result;
          rChangeCurrRing($R);
          temp->fullres = (resolvente) omAlloc0(($length + 1)*sizeof(ideal));
          for (int i = $length - 1; i >= 0; i--)
          {
             if ($r[i] != NULL)
                temp->fullres[i] = idCopy($r[i]);
          }
          temp->length = $length;
          syMinimize(temp);
          result = temp->minres;
          temp->minres = NULL;
          syKillComputation(temp, $R);
          omFreeSize((ADDRESS) temp, sizeof(ssyStrategy));
          rChangeCurrRing(origin);
          result;
       """
end

function syBetti(res::resolvente, length::Cint, R::ring)
   iv = icxx"""const ring origin = currRing;
         rChangeCurrRing($R);
         int dummy;
         intvec *iv = syBetti($res, $length, &dummy, NULL, FALSE, NULL);
         rChangeCurrRing(origin);
         return iv;
      """
   nrows = icxx"""$iv->rows();"""
   ncols = icxx"""$iv->cols();"""
   betti = icxx"""int *betti = (int *)malloc($ncols*$nrows*sizeof(int));
         for (int i = 0; i < $ncols; i++) {
            for (int j = 0; j < $nrows; j++) {
               betti[i*$nrows+j] = IMATELEM(*$iv, j+1, i+1);
            }
         }
         delete($iv);
         return &betti[0];
      """
   unsafe_wrap(Array, betti, (nrows, ncols), true)
end
