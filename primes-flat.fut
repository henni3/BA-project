-- Primes: Flat-Parallel Version
-- ==
-- compiled input { 30i64 } output { [2i64, 3i64, 5i64, 7i64, 11i64, 13i64, 17i64, 19i64, 23i64, 29i64] }
-- compiled input { 10000000i64 } auto output

-- segmented scan inclusive
let seg_scan_inc [n] 't (op: t->t->t) (ne:t) (flags: [n]bool) (imp_arr:[n]t) : [n]t =
  let (_, res) = unzip <|
    scan(\(x_flag,x)(y_flag,y) -> 
      let fl = x_flag || y_flag
      let vl = if y_flag then y else op x y
      in (fl, vl)
      ) (false, ne) (zip flags imp_arr)
  in res

-- flagarray
let mkFlagArray [m] 't (aoa_shp: [m]i64) (zero:t) (aoa_val: [m]t) : []i64 =
  let shp_rot = map(\i -> if i==0 then 0i64 else aoa_shp[i-1]) (iota m)
  let shp_scn = scan (+) 0 shp_rot
  let aoa_len = shp_scn[m-1]+aoa_shp[m-1]
  let shp_ind = map2 (\shp ind -> if shp==0 then -1 else ind) aoa_shp shp_scn
  in scatter (replicate aoa_len zero) shp_ind aoa_val

let primesFlat (n : i64) : []i64 =
  let sq_primes   = [2i64, 3i64, 5i64, 7i64]
  let len  = 8i64
  let (sq_primes, _) =
    loop (sq_primes, len) while len < n do
      -- this is "len = min n (len*len)" 
      -- but without running out of i64 bounds 
      let len = if n / len < len then n else len*len

      let mult_lens = map (\ p -> (len / p) - 1 ) sq_primes
      let flat_size = reduce (+) 0 mult_lens
      let mylen = length sq_primes 
      --------------------------------------------------------------
      -- The current iteration knowns the primes <= 'len', 
      --  based on which it will compute the primes <= 'len*len'
      -- ToDo: replace the dummy code below with the flat-parallel
      --       code that is equivalent with the nested-parallel one:
      --   let composite = map (\ p -> let mm1 = (len / p) - 1
      --                               in  map (\ j -> j * p ) (map (+2) (iota mm1))
      --                       ) sq_primes
      --   let not_primes = reduce (++) [] composite
      --
      -- Your code should compute the right `not_primes`.
      -- Please look at the lecture slides L2-Flattening.pdf to find
      --  the normalized nested-parallel version.
      -- Note that the scalar computation `mm1 = (len / p) - 1' has
      --  already been distributed and the result is stored in "mult_lens",
      --  where `p \in sq_primes`.
      -- Also note that `not_primes` has flat length equal to `flat_size`
      --  and the shape of `composite` is `mult_lens`.

      -- Flag array for iota
      let tmp = map(\i -> if i==0 then 0 else mult_lens[i-1]) (indices sq_primes)
      let shp_sc_ex = scan (+) 0 tmp
      let zeroArr = replicate flat_size 0
      let oneArr = map(\_ -> 1) shp_sc_ex-- replicate mylen 1
      let flagArr = scatter zeroArr shp_sc_ex oneArr 

      let vals = map(\f -> if f!=0 then 0 else 1) flagArr
      let iots = seg_scan_inc (+) 0 (map () flagArr) vals

      let twom = map(+2) iots
      let rpFlagArr = scatter zeroArr shp_sc_ex sq_primes
      let rp = seg_scan_inc (+) 0 flagArr rpFlagArr
      let not_primes = map(\(j,p) -> j*p) (zip twom rp)

        
     -- let not_primes = replicate flat_size 0

      -- If not_primes is correctly computed, then the remaining
      -- code is correct and will do the job of computing the prime
      -- numbers up to n!
      --------------------------------------------------------------
      --------------------------------------------------------------

       let zero_array = replicate flat_size 0i8
       let mostly_ones= map (\ x -> if x > 1 then 1i8 else 0i8) (iota (len+1))
       let prime_flags= scatter mostly_ones not_primes zero_array
       let sq_primes = filter (\i-> (i > 1i64) && (i <= n) && (prime_flags[i] > 0i8))
                              (0...len)

       in  (sq_primes, len)

  in sq_primes

-- RUN a big test with:
-- $ futhark opencl primes-flat.fut
-- $ echo "10000000" | ./primes-flat -t /dev/stderr -r 10 > /dev/null
let main (n : i64) : []i64 = primesFlat n
