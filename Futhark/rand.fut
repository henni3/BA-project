-- random module taken from Futharks webpage
module type rand ={
    type rng
    val init : i32 -> rng
    val rand : rng -> (rng, i32)
    val split : rng -> (rng,rng)
    val split_n : (n: i64) -> rng -> [n]rng
}


module lcg : rand = {
    type rng = u32
    
    def addConst : u32 = 1103515245
    def multConst : u32 = 12345
    def modConst : u32 = 1<<31

    def rand rng = 
        let rng' = (addConst * rng + multConst) % modConst
        in (rng', i32.u32 rng')
    
    
    def init (x: i32) : u32 =
        let x = u32.i32 x
        let x =((x >> 16) ^ x) * 0x45d9f3b   
        let x =((x >> 16) ^ x) * 0x45d9f3b
        let x =((x >> 16) ^ x)
        in x
    def split (rng: rng) = 
        (init (i32.u32 (rand rng).0),
         init (i32.u32 rng))

    def split_n n rng = 
        tabulate n (\i -> init (i32.u32 rng ^ i32.i64 i))
}

def rand_i32 (rng: lcg.rng) (bound: i32) =
    let (rng,x) = lcg.rand rng
    in (rng, x % bound)

let main (n:i32) : i32 = (rand_i32 (lcg.init n) 10).1
