const constants = require('./constants').constants;

const A_O_O = parseFloat(constants.A_O_O);
const p_O_O = parseFloat(constants.p_O_O);
const A_Zn_O = parseFloat(constants.A_Zn_O);
const p_Zn_O = parseFloat(constants.p_Zn_O);
const C_O_O = parseFloat(constants.C_O_O)
const k = parseFloat(constants.k);
const c = parseFloat(constants.c_force);
const m_Zn = parseFloat(constants.m_Zn);
const m_O = parseFloat(constants.m_O);
const charge_f = -2.0;
const charge_s = 2.0;

function initialization(N, coordinateOxygen, coordinateZunk) {
  var forceOX = forceZnX = new Array(N).fill(0.0);
  var forceOY = forceZnY = new Array(N).fill(0.0);
  var forceOZ = forceZnZ = new Array(N).fill(0.0);
  var accelerationOX = accelerationZnX = new Array(N).fill(0.0);
  var accelerationOY = accelerationZnY = new Array(N).fill(0.0);
  var accelerationOZ = accelerationZnZ = new Array(N).fill(0.0);
  setInterval(function () {
    force_O_O_sum(N, coordinateOxygen);
    force_O_Zn_sum(N, coordinateOxygen, coordinateZunk);
    acceleration(N);
  }, 1000);
}

function force_O_O(r) {
  force = k * (charge_f * charge_f) / (r * r) + c * (A_O_O / p_O_O) * Math.exp(-Math.abs(r) / p_O_O) + c * 6 * C_O_O / Math.pow(Math.abs(r), 7.0);
  return force;
}

function force_f_s(r) {
  force = k * (charge_s * charge_f) / (r * r) + c * (A_Zn_O / p_Zn_O) * Math.exp(-Math.abs(r) / p_Zn_O);
  return force;
}

function force_O_O_sum(N, coordinateOxygen) {
  let dx = dy = dz = r_O_O = 0;
  let f_o_o = 0;
  for (let i = 0; i < N; i++) {
    for (let j = i + 1; i < N; j++) {
      dx = coordinateOxygen.x[i] - coordinateOxygen.x[j];
      dy = coordinateOxygen.y[i] - coordinateOxygen.y[j];
      dz = coordinateOxygen.z[i] - coordinateOxygen.z[j];
      r_O_O = Math.sqrt(dx * dx + dy * dy + dz * dz);
      f_O_O = force_O_O(r_O_O);
      forceOX[i] = forceOX[i] + f_O_O * (dx / r_O_O);
      forceOX[j] = forceOX[j] - f_O_O * (dx / r_O_O);
      forceOY[i] = forceOY[i] + f_O_O * (dy / r_O_O);
      forceOY[j] = forceOY[j] - f_O_O * (dy / r_O_O);
      forceOZ[i] = forceOZ[i] + f_O_O * (dz / r_O_O);
      forceOZ[j] = forceOZ[j] - f_O_O * (dz / r_O_O);
    }
  }
}

function force_O_Zn_sum(N, coordinateOxygen,coordinateZunk) {
  let dx = dy = dz = r_Zn_O = 0;
  let f_Zn_O = 0;
  for (let i = 0; i < N; i++) {
    for (let j = 0; i < N; j++) {
      dx = coordinateOxygen.x[i] - coordinateZunk.x[j];
      dy = coordinateOxygen.y[i] - coordinateZunk.y[j];
      dz = coordinateOxygen.z[i] - coordinateZunk.z[j];
      r_Zn_O = Math.sqrt(dx * dx + dy * dy + dz * dz);
      f_Zn_O = force_O_O(r_O_O);
      forceOX[i] = forceOX[i] + f_Zn_O * (dx / r_Zn_O);
      forceZnX[j] = forceZnX[j] - f_Zn_O * (dx / r_Zn_O);
      forceOY[i] = forceOY[i] + f_Zn_O * (dy / r_Zn_O);
      forceZnY[j] = forceZnY[j] - f_Zn_O * (dy / r_Zn_O);
      forceOZ[i] = forceOZ[i] + f_Zn_O * (dz / r_Zn_O);
      forceZnZ[j] = forceZnZ[j] - f_Zn_O * (dz / r_Zn_O);
    }
  }
}

function acceleration(N){
  for(let i = 0; i < N; i++){
    accelerationOX[i] = forceOX[i]/m_O;
    accelerationOY[i] = forceOY[i]/m_O;
    accelerationOZ[i] = forceOZ[i]/m_O;
    accelerationZnX[i] = forceZnX[i]/m_Zn;
    accelerationZnY[i] = forceZnY[i]/m_Zn;
    accelerationZnZ[i] = forceZnZ[i]/m_Zn;
  }
}

function newCoordinates(N, coordinateOxygen, coordinateZunk) {
  
}

module.exports = {
  initialization: initialization,
}
