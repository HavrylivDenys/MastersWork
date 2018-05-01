const fs = require('fs');
const parseString = require('xml2js').parseString;
const util = require('util');

function readAcxFile(){	
	fs.readFile(__dirname + '/startingData/data.acx','utf8', (err, data)=>{
		parseString(data, (err, result)=> {
			parseToObject(result);
		});			
	});
}

function parseToObject(data) {
	let coord3D = {
		Number: 15,
		Oxygen: {
			x: [],
			y: [],
			z: [],
		}, 		
		Zunk: {
			x: [],
			y: [],
			z: [],
		}
	};
	for(let i = 0; i < coord3D.Number * 2; i++){
		if(data['chemistryDoc'].atom[i]['$'].atomicNumber == 8){
			coord3D.Oxygen.x.push(parseFloat(data['chemistryDoc'].atom[i]['$'].coord3DX));
			coord3D.Oxygen.y.push(parseFloat(data['chemistryDoc'].atom[i]['$'].coord3DY));
			coord3D.Oxygen.z.push(parseFloat(data['chemistryDoc'].atom[i]['$'].coord3DZ));
		}else if(data['chemistryDoc'].atom[i]['$'].atomicNumber == 30){
			coord3D.Zunk.x.push(parseFloat(data['chemistryDoc'].atom[i]['$'].coord3DX));
			coord3D.Zunk.y.push(parseFloat(data['chemistryDoc'].atom[i]['$'].coord3DY));
			coord3D.Zunk.z.push(parseFloat(data['chemistryDoc'].atom[i]['$'].coord3DZ));
		}		
	}
	writeToFile(coord3D);	
}

function writeToFile(startingData){
	fs.writeFile(__dirname + '/startingData/starting_data.json', JSON.stringify(startingData, null, '\t'), (err) =>{
		if(!err){
			console.log('All good');
		} 
	});
}

readAcxFile();

