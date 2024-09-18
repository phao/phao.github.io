function pol_degree(p) {
  return p.length - 1;
}

function pol_is0(p) {
  return p.length == 1 && p[0].isZero();;
}

function pol_ith(p, i) {
  return i < p.length ? p[i] : bigRat.zero;
}

function pol_normalize(p) {
  if (p.length == 0) {
    p.push(bigRat.zero);
    return;
  }
  
  for (;;) {
    var n = p.length;
    if (n == 1) {
      return;
    }
    if (p[n-1].isZero()) {
      p.pop();
    }
    else {
      return;
    }
  }
}

function pol_mul(p1, p2) {
  if (pol_is0(p1)) {
    return [bigRat.zero];
  }
  else if (pol_is0(p2)) {
    return [bigRat.zero];
  }
  
  var n1 = pol_degree(p1);
  var n2 = pol_degree(p2);
  var m = n1+n2;
  var pol_out = [];
  for (var i = 0; i <= m; i++) {
    var a = bigRat.zero;
    for (var j = 0; j <= i; j++) {
      if (j <= n1 && i-j <= n2) {
        a = a.add(pol_ith(p1, j).multiply(pol_ith(p2, i-j)));
      }
    }
    pol_out.push(a);
  }
  
  pol_normalize(pol_out);
  return pol_out;
}

function pol_add(p1, p2) {
  if (pol_is0(p1)) {
    return p2;
  }
  else if (pol_is0(p2)) {
    return p1;
  }
  
  var n1 = pol_degree(p1);
  var n2 = pol_degree(p2);
  var pol_out = [];
  var m = n1 > n2 ? n1 : n2;
  for (var i = 0; i <= m; i++) {
    pol_out.push(pol_ith(p1, i).add(pol_ith(p2, i)));
  }
  pol_normalize(pol_out);
  return pol_out;
}

function pol_scalar_mul(x, p) {
  if (x.isZero()) {
    return [bigRat.zero];
  }
  
  var pol_out = [];
  for (var i = 0; i < p.length; i++) {
    pol_out.push(p[i].times(x));
  }
  
  return pol_out;
}

function pol_interpolate(points) {
  var x = function (i) {
    return points[i][0];
  };
  
  var y = function (i) {
    return points[i][1];
  };
  
  var interpol = [bigRat.zero];
  
  for (var i = 0; i < points.length; i++) {
    var interpol_term = [bigRat.one];
    var scalar = bigRat.one;
    
    for (var j = 0; j < points.length; j++) {
      if (j == i) {
        continue;
      }
      scalar = scalar.multiply(x(i).minus(x(j)));
      interpol_term = pol_mul(interpol_term,
                              [x(j).times(bigRat.minusOne), 1]);
    }
     
    interpol_term = pol_scalar_mul(y(i).divide(scalar), interpol_term);
    interpol = pol_add(interpol, interpol_term);
  }
  
  pol_normalize(interpol);
  return interpol;
}

function sum_formula(power) {
  var points = [];
  var sum_so_far = bigRat.zero;
  for (var i = 1; i <= power+2; i++) {
    var n = bigRat(i);
    sum_so_far = sum_so_far.add(n.pow(power));
    points.push([n, sum_so_far]);
  }
  return pol_interpolate(points);
}

function rat_latex(r) {
  out = "";
  out = "{";
  out += r.num.toString();
  den = r.denom;
  if (den.compare(1) != 0) {
    out += " \\over ";
    out += den.toString();
  }
  out += "}";
  return out;
}

function pol_latex(pol, variable) {
  out = "";
  
  if (pol_is0(pol)) {
    return "0";
  }
  
  put_plus = true;
  
  if (pol[0].compare(0) != 0) {
    out += rat_latex(pol[0]);
  }
  else {
    put_plus = false;
  }
  
  for (var i = 1; i < pol.length; i++) {
    if (pol[i].compare(0) == 0) {
      continue;
    }
    
    if (put_plus) {
      if (pol[i].compare(0) < 0) {
        out += " - ";
      }
      else {
        out += " + ";
      }
    }
    
    var abs_coef = pol[i].abs()
    if (abs_coef.compare(1) != 0) {
      out += rat_latex(abs_coef);
    }
    out += "{" + variable + "}";
    if (i > 1) {
      out += "^{" + i.toString() + "}";
    }
    put_plus = true;
  }
  
  return out;
}

function pol_sub(p1, p2) {
  return pol_add(p1, pol_scalar_mul(bigRat.minusOne, p2));
}

function sum_formula_rec(power) {
  BR0 = bigRat.zero;
  BR1 = bigRat.one;
  BRm1 = bigRat.minusOne;

  var binomial_coef_data = Array(power+2)
  for (var i = 0; i <= power+1; i++) {
    binomial_coef_data[i] = Array(i+1)
    binomial_coef_data[i][0] = BR1;
    binomial_coef_data[i][i] = BR1;
  }

  var binomial_coef = function (i,j) {
    if (binomial_coef_data[i][j] == undefined) {
      b_im1_jm1 = binomial_coef(i-1,j-1);
      b_im1_j = binomial_coef(i-1,j);
      binomial_coef_data[i][j] = b_im1_jm1.add(b_im1_j);
    }
    return binomial_coef_data[i][j];
  }

  np1_to_k = Array(power+2);
  np1_to_k[0] = [BR1];
  np1 = [BR1, BR1];
  for (var k = 1; k <= power+1; k++) {
    np1_to_k[k] = pol_mul(np1_to_k[k-1], np1);
  }

  pol_1 = [BR1];
  the_pols = Array(power+1);
  the_pols[0] = [BR0, BR1];

  for (var k = 1; k <= power; k++) {
    p = pol_sub(np1_to_k[k+1], pol_1);
    for (var j = 2; j <= k+1; j++) {
      p = pol_sub(p, pol_scalar_mul(binomial_coef(k+1, j), the_pols[k+1-j]));
    }
    inv_kp1_choose_1 = BR1.divide(binomial_coef(k+1, 1));
    the_pols[k] = pol_scalar_mul(inv_kp1_choose_1, p);
  }

  return the_pols[power];
}
