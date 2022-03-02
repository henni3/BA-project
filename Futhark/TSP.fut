-- mkFlagArray is taken from PMPH lecture notes p. 48
let mkFlagArray 't [m] 
        (aoa_shp: [m]i32) (zero: t)
        (aoa_val: [m]t) : []t =
    let shp_rot = map (\i -> if i == 0 then 0
                             else aoa_shp[i-1]
                      ) (iota m)
    let shp_scn = scan (+) 0 shp_rot    
    let aoa_len = shp_scn[m-1]+ aoa_shp[m-1] |> i64.i32
    let shp_ind = map2 (\shp ind -> 
                        if shp == 0 then -1
                        else ind
                        ) aoa_shp shp_scn
    in scatter (replicate aoa_len zero) (map i64.i32 shp_ind) aoa_val

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
let findMinChange [m] [n] (dist : [m]i32) (tour : [n]i32) (cities : i32) : (i32, i32, i32) =
    let totIter = ((cities-1)*(cities-2))/2 |> i64.i32
    let len = i64.i32 (cities-2)
    let aoa_val = replicate len 1i32
    let temp = map (+1) (iota len) |> map i32.i64
    let aoa_shp = map (\i ->
                       temp[len - i - 1] 
                      ) (iota len)
    let flagArr = mkFlagArray aoa_shp 0i32 aoa_val
    let Iarr = scan (+) 0i32 flagArr |> map (\x -> x-1)
    let Jarr = segmented_scan (+) 0i32 (map bool.i32 (flagArr :> [totIter]i32)) (replicate totIter 1i32) |> map (\x -> x-1)
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
let swap [m] (i : i32) (j : i32) (tour : [m]i32) : [m]i32 =
    let minI = i+1 in 
    map i32.i64 (iota m) |>
    map(\ind ->
        if ind < minI || ind > j then
            tour[ind]
        else
            tour[j - (ind - minI)]
    ) 

let cost [n] [m] (tour : [n]i32) (distM : [m]i32) : i32 = 
        map (\i ->  distM[tour[i] * i32.i64(n-1) + tour[i+1]]
            ) (iota (n-1)) |> reduce (+) 0

let twoOptAlg [m] [n] (distM : [m]i32) (tour : [n]i32) (cities : i32) : (i32, []i32) =
    --let retTour = tour
    let init = findMinChange distM tour cities
    --let xs = findMinChange distM tour cities
    let rs = loop (i, xs, cond) = (0, tour, init)  while cond.0 < 0 do 
        (i+cond.0, swap cond.1 cond.2 xs, findMinChange distM xs cities)
    in (rs.0, rs.1)

let main : (i32, i32, []i32, i32) =
    let cities = 5
    let tour = [4,2,3,1,0,4]
    let distM = [0,4,6,8,3,
                 4,0,4,5,2,
                 6,4,0,2,3,
                 8,5,2,0,4,
                 3,2,3,4,0]
    let oldCost = cost tour distM            
    --let minTour = twoOptAlg distM tour cities
    let firstChange = findMinChange distM tour cities
    let firstTour = swap firstChange.1 firstChange.2 tour
    let newCost = cost firstTour distM
    in (oldCost, newCost, firstTour, firstChange.0)
    