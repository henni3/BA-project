-- mkFlagArray is taken from PMPH lecture notes p. 48
let mkFlagArray 't [m] 
        (aoa_shp: [m]i64) (zero: t)
        (aoa_val: [m]t) : []i64 =
    let shp_rot = map (\i -> if i == 0 then 0
                             else aoa_shp[i-1]
                      ) (iota m)
    let shp_scn = scan (+) 0 shp_rot    
    let aoa_len = shp_scn[m-1]+ aoa_shp[m-1]
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
     
-- Comparator function used in forLoops in the reduce part
let changeComparator (t1 : (i32, i32, i32)) (t2: (i32, i32, i32)) : (i32, i32, i32) =
    if t1.0 < t2.0 then t1
    else 
        if t1.0 == t2.0 then 
            if t1.1 < t2.1 then t1
            else 
                if t1.1 == t2.1 then
                    if t1.2 < t2.2 then t1
                    else t2
                else t2
        else t2

-- findMinChange is the parallel implementation of the two for loops
-- in the 2-opt move algorithm
let findMinChange [m] [n] (distM : [m]i32) (tour : [n]i32) (cities : i32) : (i32, i32, i32) =
    let totIter = ((cities-1)*(cities-2))/2
    let len = cities-2
    let aoa_val = replicate len 1i32
    let temp = map (+1) (iota len)
    let aoa_shp = map (\i ->
                       temp[len - i - 1] 
                      ) (iota len)
    let flagArr = mkFlagArray aoa_shp 0i64 aoa_val
    let Iarr = scan (+) 0i64 flagArr |> map (\x -> x-1)
    let Jarr = segmented_scan (+) 0i64 flagArr (replicate totIter 1) |> map (\x -> x-1)
    let changeArr = map (\ind -> 
                        let i = Iarr[ind]
                        let iCity = tour[i]
                        let iCityp1 = tour[i+1]
                        let j = Jarr[ind] + i + 2
                        let jCity = tour[j]
                        let jCityp1 = tour[j+1]
                        in  ((dist[iCity * cities + jCity] + 
                            dist[iCityp1 * cities + jCityp1] - 
                            (dist[iCity * cities + iCityp1] + 
                            dist[jCity * cities + jCityp1])), i, j) 
                        ) (iota totIter)
    in reduce changeComparator (2147483647, -1, -1) changeArr

-- This function swaps the two edges that produce the lowest cost
let swap [m] (i : i32) (j : i32) (cities : i32) (tour : [m]i32) : [m]i3 =
    map(\ind ->
        if ind < i || ind > j then
            tour[ind]
        else
            tour[j - (ind - i)]
    ) (iota cities + 1) 


let main (cities : i32) : []i32 =
    let cities = 5
    let distM = [0,4,6,8,3,
                 4,0,4,5,2,
                 6,4,0,2,3,
                 8,5,2,0,4,
                 3,2,3,4,0]
    
    let minChange = findMinChange distM tour cities
    in swap minChange.1 minChange.2 tour


    

    
    
