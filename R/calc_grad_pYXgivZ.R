calc_grad_pYXgivZ = function(x, y, z, object) {
  calc_grad_pYgivXZ(object = object, y = y, x = x, z = z) *
    calc_grad_pXgivZ(object = object, x = x, z = z)
}
