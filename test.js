// Molecular dynamics version 1.0.0
const fs = require('fs');
const path = require('path');
const starting_data = require('./starting_data.json');

const constants = require('./constants').constants;

let A_O_O;
let p_O_O;
let A_Zn_O;
let p_Zn_O;
let C_O_O;
let k = parseFloat(constants.k);
let c = parseFloat(constants.c_force);
let m_Zn = parseFloat(constants.m_Zn);
let m_O = parseFloat(constants.m_O);			

function setConstants(name){
	switch(name){
		case 'Berhman': 
			A_O_O = parseFloat(constants.A_O_O);
			p_O_O = parseFloat(constants.p_O_O);
			A_Zn_O = parseFloat(constants.A_Zn_O);
			p_Zn_O = parseFloat(constants.p_Zn_O);
			C_O_O = parseFloat(constants.C_O_O);
			break;
		case 'Nyberg': 
			A_O_O = parseFloat(constants.A_O_O_Nyberg);
			p_O_O = parseFloat(constants.p_O_O_Nyberg);
			A_Zn_O = parseFloat(constants.A_Zn_O_Nyberg);
			p_Zn_O = parseFloat(constants.p_Zn_O_Nyberg);
			C_O_O = parseFloat(constants.C_O_O_Nyberg);
			break;
		case 'test':
			k = 1;
			c = 1;
			m_Zn = 60910.39E6;
			m_O = 14903.34E6;
			A_O_O = parseFloat(constants.A_O_O);
			p_O_O = parseFloat(constants.p_O_O);
			A_Zn_O = parseFloat(constants.A_Zn_O);
			p_Zn_O = parseFloat(constants.p_Zn_O);
			C_O_O = parseFloat(constants.C_O_O);
			break;			
	}	
}

const charge_f = -2.0;
const charge_s = 2.0;

const assemblyName = '(ZnO)\u2081\u2085';

const O_coordinates = starting_data.Oxygen;
const Zn_coordinates = starting_data.Zunk;
const N = starting_data.Number;
const t = 1e-15;

class Atom{
	constructor(x = 0, y = 0, z = 0){
		this.x = x;
		this.y = y;
		this.z = z;
		this.x0 = x;
		this.y0 = y;
		this.z0 = z;
		this.forceX = 0;
		this.forceY = 0;
		this.forceZ = 0;
		this.accelerationX = 0.0;
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
	Oxygen[i] = new Atom(O_coordinates.x[i],O_coordinates.y[i],O_coordinates.z[i]);
}
const Zunk = new Array(N);
for(let i = 0; i < Zunk.length; i++){
	Zunk[i] = new Atom(Zn_coordinates.x[i],Zn_coordinates.y[i],Zn_coordinates.z[i]);
}

function countForce(Oxygen, Zunk) {
    let dx = 0;
    let dy = 0;
    let dz = 0;
    let dr = 0;
    for (let i = 0; i < N; i++) {
      Oxygen[i].forceX = 0;
      Oxygen[i].forceY = 0;
      Oxygen[i].forceZ = 0;
      Zunk[i].forceX = 0;
      Zunk[i].forceX = 0;
      Zunk[i].forceX = 0;
    }
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

function verle(Oxygen, Zunk){
	let a = 0; let b = 0; let c = 0;
	let d = 0; let f = 0; let e = 0;
	for (let i = 0; i < N; i++) {
		a =  Oxygen[i].x;
		b =  Oxygen[i].y;
		c =  Oxygen[i].z;
		d =  Zunk[i].x;
		f =  Zunk[i].y;
		e =  Zunk[i].z; 	 		
		Oxygen[i].x = 2 * Oxygen[i].x - Oxygen[i].x0 + Oxygen[i].accelerationX * t * t;
		Oxygen[i].y = 2 * Oxygen[i].y - Oxygen[i].y0 + Oxygen[i].accelerationY * t * t;
		Oxygen[i].z = 2 * Oxygen[i].z - Oxygen[i].z0 + Oxygen[i].accelerationZ * t * t;
		Zunk[i].x = 2 * Zunk[i].x - Zunk[i].x0 + Zunk[i].accelerationX * t * t;
		Zunk[i].y = 2 * Zunk[i].y - Zunk[i].y0 + Zunk[i].accelerationY * t * t;
		Zunk[i].z = 2 * Zunk[i].z - Zunk[i].z0 + Zunk[i].accelerationZ * t * t;
		Oxygen[i].x0 = a;
		Oxygen[i].y0 = b;
		Oxygen[i].z0 = c;
		Zunk[i].x0 = d;
		Zunk[i].y0 = f;
		Zunk[i].z0 = e;
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
function workFlow(name){
	let cycle = 0;
	setConstants(name);
	for(let i = 0; i <= 1e7; i++){
		countForce(Oxygen, Zunk);
		countAcceleration(Oxygen, Zunk);
		if(i == 0){
			countSpeed(Oxygen, Zunk);
			countCoordinates(Oxygen, Zunk);
		}else {
			verle(Oxygen,Zunk);
		}		
		if((i == 0 || i == 1e4 || i == 1e5 || i == 1e7)){
			writeToAFile(name, cycle);
			cycle++;
		}		
	}	
}

function writeToAFile(name, cycle){
	let dir = {dir: __dirname, base: name};    
    let dirName = path.format(dir);
    let file = {dir: dirName, base: `${cycle}.xyz`};
    let fileName = path.format(file);
	if(!fs.existsSync(dirName)){
		fs.mkdirSync(dirName);
	}
	try{
		let O;
		let Zn;				
		fs.appendFileSync(fileName, `${N}\n`);
		fs.appendFileSync(fileName, `${assemblyName}\n`);
		for(let j = 0; j < N; j++){
			O = `O`.padEnd(10,' ') + `${Oxygen[j].x.toFixed(9)}`.padEnd(15, ' ') 
			+ `${Oxygen[j].y.toFixed(9)}`.padEnd(15, ' ') + `${Oxygen[j].z.toFixed(9)}`.padEnd(15, ' ') + '\n';
			Zn = `Zn`.padEnd(10,' ') + `${Zunk[j].x.toFixed(9)}`.padEnd(15, ' ') 
			+ `${Zunk[j].y.toFixed(9)}`.padEnd(15, ' ') + `${Zunk[j].z.toFixed(9)}`.padEnd(15, ' ') + '\n';
			fs.appendFileSync(fileName, O);
			fs.appendFileSync(fileName, Zn);
		}		
	}catch(err){
		console.log(err);			
	}
}

function start(){
	workFlow('Berhman');
	workFlow('Nyberg');
	workFlow('test');
}

// Test functions

function countTime(){
	let start,end;	
	let O, Z;
	let time = 0;
	// lets make a copy of objects to make counting clear function	
	O = JSON.parse(JSON.stringify(Oxygen));
	Z = JSON.parse(JSON.stringify(Zunk));
	start = new Date();
	for(let i = 0; i < 100; i++){			
		countForce(O, Z);
		countAcceleration(O, Z);
		countSpeed(O, Z);
		countCoordinates(O, Z);	
	}
	end = new Date();
	time = (end - start)/100;
	return time;
} // counting time for cicle
const startO = JSON.parse(JSON.stringify(Oxygen));
function print(){
	for(let i = 0; i < N; i++){
		console.log(startO[i].x - Oxygen[i].x, i);
	}
} // test function


function test(){
	setConstants('Berhman');
	let cycle = 0;
	for(let i = 0; i <= 1e4; i++){
		countForce(Oxygen, Zunk);
		countAcceleration(Oxygen, Zunk);
		if(i == 0){
			countSpeed(Oxygen, Zunk);
			countCoordinates(Oxygen, Zunk);
		}else {
			verle(Oxygen,Zunk);
		}		
		if((i == 0 || i == 1e4 || i == 1e5 || i == 1e7)){
			writeToAFile('Berhman', cycle);			
			cycle++;
		}		
	}	
}

test();
print();