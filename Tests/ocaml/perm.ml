
let ctz x =
  assert (x <> 0);
  let x = ref x in
  let n = ref 62 in
  let y = !x lsl 32 in if y <> 0 then (n := !n - 32 ; x := y);
  let y = !x lsl 16 in if y <> 0 then (n := !n - 16 ; x := y);
  let y = !x lsl  8 in if y <> 0 then (n := !n -  8 ; x := y);
  let y = !x lsl  4 in if y <> 0 then (n := !n -  4 ; x := y);
  let y = !x lsl  2 in if y <> 0 then (n := !n -  2 ; x := y);
  let y = !x lsl  1 in if y <> 0 then (n := !n -  1 ; x := y);
  !n
;;

let compute n m = 
  let k = ref 0 in
  let u = ref ((1 lsl n) - 1) in
  let v = Array.make 10000 0 in
  while !u < (1 lsl (n+m)) do
    v.(!k) <- !u;
    incr k;
    let t = !u lor (!u-1) in
    let t' = t+1 in
    let t'' = ( (((lnot t) land t')-1) lsr (ctz !u + 1) ) in
    u := t' lor t'';
  done;
  v
;;

let bin_of_int d =
  if d < 0 then invalid_arg "bin_of_int" else
  if d = 0 then "0" else
  let rec aux acc d =
    if d = 0 then acc else
    aux (string_of_int (d land 1) :: acc) (d lsr 1)
  in
  String.concat "" (aux [] d)


let () =
  compute 4 4
  |> Array.iteri (fun i v -> Printf.printf "%d %s\n" i (bin_of_int v)) 


