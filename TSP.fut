-- mkFlagArray is taken from PMPH lecture notes p. 48
let mkFlagArray 't [m] 
        (aoa_shp: [m]i32) (zero: t)
        (aoa_val: [m]t) : []i32 =
    let shp_rot = map (\i -> if i == 0 then 0
                             else aoa_shp[i-1]
                      ) (iota m)
    let shp_scn = scan (+) 0 shp_rot
    let aoa_len = shp_scn[m-1]+aoa_shp[m-1]
    let shp_ind = map2 (\shp ind -> 
                        if shp == 0 then -1
                        else ind
                        ) aoa_shp shp_scn
    in scatter (replicate aoa_len zero) shp_ind aoa_val

-- segmented_scan is taken from PMPH Futhark code 
let segmented_scan [n] 't (op: t -> t -> t) (ne: t)
                          (flags: [n]bool) (arr: [n]t) : [n]t =
    let (_, res) = unzip <|
        scan (\(x_flag,x) (y_flag,y) ->
                let fl = x_flag || y_flag
                let vl = if y_flag then y else op x y
                in  (fl, vl)
            ) (false, ne) (zip flags arr)
    in  res
     



let main (n : i32) : []i32
    let totIter = ((n-1)*(n-2))/2
    let len = n-2
    let aoa_val = replicate len 1
    let temp = map (+1) (iota len)
    let aoa_shp = map (\i ->
                       temp[len - i - 1] 
                      ) (iota len)
    let flagArr = mkFlagArray aoa_shp 0 aoa_val
    let Iarr = scan (+) 0 flagArr |> map (\x -> x-1)
    let Jarr = segmented_scan (+) 0 flagArr (replicate totIter 1) |> map (\x -> x-1)
    Iarr