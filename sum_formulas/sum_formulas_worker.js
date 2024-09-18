importScripts("BigInt_BigRat.min.js", "sum_formulas.js");

onmessage = function (msg) {
  p = msg.data;
  postMessage(pol_latex(sum_formula_rec(p), "n"));
}
