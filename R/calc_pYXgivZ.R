# Calculate probabilities/densities from model of Y,X|Z
calc_pYXgivZ = function(object, y, x, z) {
  calc_pYgivXZ(object = object, y = y, x = x, z = z) *
    calc_pXgivZ(object = object, x = x, z = z)
}
