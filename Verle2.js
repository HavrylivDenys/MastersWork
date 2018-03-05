const constants = require('./constants').constants;

const A_O_O = parseFloat(constants.A_O_O);
const p_O_O = parseFloat(constants.p_O_O);
const A_Zn_O = parseFloat(constants.A_Zn_O);
const p_Zn_O = parseFloat(constants.p_Zn_O);
const C_O_O = parseFloat(constants.C_O_O)
const k = parseFloat(constants.k);
const c = parseFloat(constants.c_force);
const charge_O = -2.0;
const charge_Zn = 2.0;

function force_O_O(r) {
  force = k * (charge_O * charge_O) / (r * r) + c * (A_O_O / p_O_O) * Math.exp(-Math.abs(r) / p_O_O) + c * 6 * C_O_O / Math.pow(Math.abs(r), 7.0);
  return force;
}

function force_Zn_O(r) {
  force = k * (charge_Zn * charge_O) / (r * r) + c * (A_Zn_O / p_Zn_O) * Math.exp(-Math.abs(r) / p_Zn_O);
  return force;
}

function verleMethod(first_particle, second_particle, N) {
  let dx_f_f = dy_f_f = dz_f_f = r_f_f = 0.0;
  let dx_s_s = dy_s_s = dz_s_s = r_f_s = 0.0;
  let force_f_f_x = force_f_s_x = new Array(4);
  let force_f_f_y = force_f_s_y = new Array(4);
  let force_f_f_z = force_f_s_z = new Array(4);
  let f_o = f_zn = 0.0;

  // For O-O interaction
  for (let i = 0; i < N; i++) {
    for (let j = i + 1; i < N; j++) {
      dx_f_f = first_particle.x[i] - first_particle.x[j];
      dy_f_f = first_particle.y[i] - first_particle.y[j];
      dz_f_f = first_particle.z[i] - first_particle.z[j];

      r_f_f = Math.sqrt(dx_f_f * dx_f_f + dy_f_f * dy_f_f + dz_f_f * dz_f_f);
      f_o = force_O_O(r_f_f);

      force_f_f_x[i] = force_f_f_x[i] + f_o * (dx_f_f / r_f_f);
      force_f_f_x[j] = force_f_f_x[j] - f_o * (dx_f_f / r_f_f);
      force_f_f_y[i] = force_f_f_y[i] + f_o * (dy_f_f / r_f_f);
      force_f_f_y[j] = force_f_f_y[j] - f_o * (dy_f_f / r_f_f);
      force_f_f_z[i] = force_f_f_z[i] + f_o * (dz_O_O / r_f_f);
      force_f_f_z[j] = force_f_f_z[j] - f_o * (dz_O_O / r_f_f);
    }
  }
  // For Zn-O interaction
  for (let i = 0; i < N; i++) {
    for (let j = 0; i < N; j++) {
      dx_s_s = first_particle.x[i] - second_particle.x[j];
      dy_s_s = first_particle.y[i] - second_particle.y[j];
      dz_s_s = first_particle.z[i] - second_particle.z[j];

      r_f_s = Math.sqrt(dx_s_s * dx_s_s + dy_s_s * dy_s_s + dz_s_s * dz_s_s);
      f_zn = force_Zn_O(r_f_s);

      force_f_s_x[i] = force_f_s_x[i] + f_o * (dx_s_s / r_f_s);
      force_f_s_x[j] = force_f_s_x[j] - f_o * (dx_s_s / r_f_s);
      force_f_s_y[i] = force_f_s_y[i] + f_o * (dy_s_s / r_f_s);
      force_f_s_y[j] = force_f_s_y[j] - f_o * (dy_s_s / r_f_s);
      force_f_s_z[i] = force_f_s_z[i] + f_o * (dz_s_s / r_f_s);
      force_f_s_z[j] = force_f_s_z[j] - f_o * (dz_s_s / r_f_s);
    }
  }
}

module.exports = {
  verleMethod: verleMethod
}
