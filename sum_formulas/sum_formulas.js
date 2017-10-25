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
    out += rat_latex(rat_latex(pol[0]));
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
