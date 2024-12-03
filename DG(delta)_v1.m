(* ::Package:: *)

(* ::Input::Initialization:: *)
path = "C:/Users/User/Desktop/master/RL/CFT_bootstrap/LP_newEnv/LP_data/isingDGn_a.csv";
(*\:8bfb\:53d6CSV\:6587\:4ef6*)
data=Import[path,"CSV"];
processedData = data/. s_String:>ToExpression[s];
(*Call constraintnum by their derivative order
f: first, l: last*)
DG[delta_, fconstraintnum_, lconstraintnum_]:=Table[N[processedData[[i]][[2]]/.{d->delta, z->1/2}, 100]/(Factorial[processedData[[i]][[1]]]/.{d->delta, z->1/2}), {i,  (fconstraintnum+1)/2+1,(lconstraintnum+1)/2+1}];
