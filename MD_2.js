// Molecular dynamics version 1.0.0
const fs = require('fs');
const starting_data = require('./starting_data.json');

const constants = require('./const_v_2').constants;

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
const t = 1;

class Molecule{
	constructor(x = 0, y = 0, z = 0){
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
		// Counting force between O-O molecules
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
		// Counting force between Zn-O molecules
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
function countAcceleration(Oxygen, Zunk) {
	for (let i = 0; i < N; i++) {
			Oxygen[i].accelerationX = Oxygen[i].forceX / m_O;
			Oxygen[i].accelerationY = Oxygen[i].forceY / m_O;
			Oxygen[i].accelerationZ = Oxygen[i].forceZ / m_O;
			Zunk[i].accelerationX = Zunk[i].forceX / m_Zn;
			Zunk[i].accelerationY = Zunk[i].forceY / m_Zn;
			Zunk[i].accelerationZ = Zunk[i].forceZ / m_Zn;					
	}
}
function countSpeed(Oxygen, Zunk) {
	for (let i = 0; i < N; i++) {
			Oxygen[i].speedX = Oxygen[i].speedX + Oxygen[i].accelerationX * t;
			Oxygen[i].speedY = Oxygen[i].speedY + Oxygen[i].accelerationY * t;
			Oxygen[i].speedZ = Oxygen[i].speedY + Oxygen[i].accelerationY * t;
			Zunk[i].speedX = Zunk[i].speedX + Zunk[i].accelerationX * t;
			Zunk[i].speedY = Zunk[i].speedY + Zunk[i].accelerationY * t;
			Zunk[i].speedZ = Zunk[i].speedY + Zunk[i].accelerationY * t;				
	}
}
function countCoordinates(Oxygen, Zunk) {
	for (let i = 0; i < N; i++) {
			Oxygen[i].x = Oxygen[i].x + Oxygen[i].speedX * t;
			Oxygen[i].y = Oxygen[i].y + Oxygen[i].speedY * t;
			Oxygen[i].z = Oxygen[i].z + Oxygen[i].speedY * t;
			Zunk[i].x = Zunk[i].x + Zunk[i].speedX * t;
			Zunk[i].y = Zunk[i].y + Zunk[i].speedY * t;
			Zunk[i].z = Zunk[i].z + Zunk[i].speedY * t;				
	}
}
function force_O_O(r) {
  force = k * (charge_f * charge_f) / (r * r) + c * (A_O_O / p_O_O) * Math.exp(-Math.abs(r) / p_O_O) - c * 6 * C_O_O / Math.pow(Math.abs(r), 7.0);
  return force;
}
function force_Zn_O(r) {
  force = k * (charge_s * charge_f) / (r * r) + c * (A_Zn_O / p_Zn_O) * Math.exp(-Math.abs(r) / p_Zn_O);
  return force;
}

// const interval = setInterval(() =>{	
// 	countForce(Oxygen, Zunk);
// 	countAcceleration(Oxygen, Zunk);
// 	countSpeed(Oxygen, Zunk);
// 	countCoordinates(Oxygen, Zunk);	
// }, 10); // real time of execution on my machine is 6-7 ms.

// setTimeout(() =>{
// 	clearInterval(interval);
// 	// print();
// }, 200000);
function workFlow(){
	for(let i = 0; i <= 1e5; i++){
		countForce(Oxygen, Zunk);
		countAcceleration(Oxygen, Zunk);
		countSpeed(Oxygen, Zunk);
		countCoordinates(Oxygen, Zunk);
		if(!(i % 100000)){
			try{				
				for(let j = 0; j < N; j++){			
					fs.appendFileSync('results2.txt', `Oxygen x: ${Oxygen[j].x.toFixed(9)} y: ${Oxygen[j].y.toFixed(9)} z: ${Oxygen[j].z.toFixed(9)}\nZunk x: ${Zunk[j].x.toFixed(9)} y: ${Zunk[j].y.toFixed(9)} z: ${Zunk[j].z.toFixed(9)} cycle: ${i}\n`,);
				}
			}catch(err){
				console.log(err);			
			}finally{
				fs.appendFileSync('results2.txt', '\n');		
			}			
		}
	}	
}

function countTime(){
	let start,end;	
	let O, Z;
	let time = 0;
	// lets make a copy of objects to make counting clear function	
	O = JSON.parse(JSON.stringify(Oxygen));
	Z = JSON.parse(JSON.stringify(Zunk));
	for(let i = 0; i < 100; i++){
		start = new Date();	
		countForce(O, Z);
		countAcceleration(O, Z);
		countSpeed(O, Z);
		countCoordinates(O, Z);
		end = new Date();
		if(time < (end - start)){
			time = end - start;
		}
	}	
	return time;
} // counting time for cicle

const startO = JSON.parse(JSON.stringify(Oxygen));
function print(){
	for(let i = 0; i < N; i++){
		console.log(startO[i].x - Oxygen[i].x);		
	}
} // test function

workFlow();
print();