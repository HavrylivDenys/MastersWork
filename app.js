const MolecularDynamics = require('./MolecularDynamics');

const starting_data = require('./starting_data.json');

const O_coordinates = starting_data.Oxygen;
const Zn_coordinates = starting_data.Zunk;
const N = starting_data.Number;

MolecularDynamics.initialization(N,O_coordinates,Zn_coordinates);
