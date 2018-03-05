// Molecular dynamics version 1.0.0

const starting_data = require('./starting_data.json');

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

const O_coordinates = starting_data.Oxygen;
const Zn_coordinates = starting_data.Zunk;
const N = starting_data.Number;

class Molecule{
	constructor(x = 0,y = 0,z = 0){
		this.x = x;
		this.y = y;
		this.z = z;
		this.forceX = 0;
		this.forceY = 0;
		this.forceZ = 0;
		this.accelerationX = 0;
		this.accelerationY = 0;
		this.accelerationZ = 0;
		this.speedX = 0;
		this.speedY = 0;
		this.speedZ = 0;
	}
	setCoordinates(x,y,z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	print(){
		console.log(this.x, this.y, this.z);
	}
}

const Oxygen = new Array(N);
for(let i = 0; i < Oxygen.length; i++){
	Oxygen[i] = new Molecule(O_coordinates.x[i],O_coordinates.y[i],O_coordinates.z[i]);
}
const Zunk = new Array(N);
for(let i = 0; i < Zunk.length; i++){
	Zunk[i] = new Molecule(Zn_coordinates.x[i],Zn_coordinates.y[i],Zn_coordinates.z[i]);
}
function countForce(Oxygen, Zunk) {
    let dx = dy = dz = dr = 0;
    for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
        	dx = Oxygen[i].x - Oxygen[j].x;
        	dy = Oxygen[i].y - Oxygen[j].y;
        	dz = Oxygen[i].z - Oxygen[j].z;
        	dr = Math.sqrt(dx * dx + dy * dy + dz * dz);
        	Oxygen[i].forceX = Oxygen[i].forceX + (dx/dr) * force_O_O(dr);        	
        	Oxygen[j].forceX = Oxygen[j].forceX - (dx/dr) * force_O_O(dr);
        	Oxygen[i].forceY = Oxygen[i].forceY + (dy/dr) * force_O_O(dr);
        	Oxygen[j].forceY = Oxygen[j].forceY - (dy/dr) * force_O_O(dr);
        	Oxygen[i].forceZ = Oxygen[i].forceZ + (dz/dr) * force_O_O(dr);
        	Oxygen[j].forceZ = Oxygen[j].forceZ - (dz/dr) * force_O_O(dr);
        }
    }
    dx = dy = dz = dr = 0;

    for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
        	dx = Oxygen[i].x - Zunk[j].x;
        	dy = Oxygen[i].y - Zunk[j].y;
        	dz = Oxygen[i].z - Zunk[j].z;
        	dr = Math.sqrt(dx * dx + dy * dy + dz * dz);
        	Oxygen[i].forceX = Oxygen[i].forceX + (dx/dr) * force_Zn_O(dr);        	
        	Zunk[j].forceX = Zunk[j].forceX - (dx/dr) * force_Zn_O(dr);
        	Oxygen[i].forceY = Oxygen[i].forceY + (dy/dr) * force_Zn_O(dr);
        	Zunk[j].forceY = Zunk[j].forceY - (dy/dr) * force_Zn_O(dr);
        	Oxygen[i].forceZ = Oxygen[i].forceZ + (dz/dr) * force_Zn_O(dr);
        	Zunk[j].forceZ = Zunk[j].forceZ - (dz/dr) * force_Zn_O(dr);
        }
    }
}

function force_O_O(r) {
  force = k * (charge_f * charge_f) / (r * r) + c * (A_O_O / p_O_O) * Math.exp(-Math.abs(r) / p_O_O) + c * 6 * C_O_O / Math.pow(Math.abs(r), 7.0);
  return force;
}

function force_Zn_O(r) {
  force = k * (charge_s * charge_f) / (r * r) + c * (A_Zn_O / p_Zn_O) * Math.exp(-Math.abs(r) / p_Zn_O);
  return force;
}
countForce(Oxygen, Zunk);

console.log(Zunk);