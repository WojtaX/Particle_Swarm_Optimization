double Uniform_Mersenne(double Low, double Up) {
	random_device seed;
	mt19937 gen(seed());
	uniform_real_distribution<double> dist(Low, Up);
	return dist(gen);
}

void Swarm::Sorting() {//Sortowanie cz¹stek ze wzglêdu na pozycjê. 
	for (int j = 0; j < Particles.size(); j++) {
		for (int i = 0; i < Particles.size() - 1; i++) {
			if (Particles[i].Value > Particles[i + 1].Value)
				swap(Particles[i], Particles[i + 1]);

		}
	}

}

double Swarm::FunctionCalculating(vector<double>x) {          //Funkcja do optymalizacji 
	double ValueForParticle = 0;
	for (int i = 0; i < x.size(); i++) {
		//ValueForParticle += cos(Position[i])*cos(3*Position[i])*cos(5*Position[i]); 
		//ValueForParticle += (1 / (0.1 + abs(x[i]))*(1 + cos(x[i])));
		// ValueForParticle +=(-1*(Position[i] - 10)*(Position[i] - 10) + 100);
		ValueForParticle += (10 * (sin(0.2*(x[0] - 10)) + 2) * (sin(0.5*(x[i] - 10))) * (sin(0.5*(x[i] - 10))) / abs(x[i] - 10)) + 1;
	} 
	return ValueForParticle;
}
double Swarm::FunctionCalculating(double x) {          //Funkcja do optymalizacji 
	double ValueForParticle = 0;
	for (int i = 0; i < Dimension; i++) {
		//ValueForParticle += cos(Position)*cos(3*Position)*cos(5*Position); 
		//ValueForParticle += (1 / (0.1 + abs(x))*(1 + cos(x)));
		//ValueForParticle += (-1 * (Position - 10)*(Position - 10) + 1240);
		ValueForParticle += (10 * (sin(0.2*(x - 10)) + 2) * (sin(0.5*(x - 10))) * (sin(0.5*(x - 10))) / abs(x - 10)) + 1;
		
	}
	return ValueForParticle;
}

Particle::Particle(int Dimension) {   //konstruktor cz¹stki

	for (int i = 0; i < Dimension; i++) {
		double NewRandomParticle = Uniform_Mersenne(LOW, UP); //generowanie nowej cz¹stki 
		Position.push_back(NewRandomParticle);
		Velocity.push_back(0);
		ParticleBest.push_back(0);
	}
}

void  Swarm::SwarmGenerating() {  //Tworzenie roju o rozmiarze "SwarmSize", przy uzyciu wartosci losowanych Uniform_Mersenne 
									//z przedzialu [LOW , UP]
	for (int i = 0; i < SwarmSize; i++) {
		Particle Wojtek(Dimension); //tworzenie nowej cz¹stki
		Particles.push_back(Wojtek);
	}
}

void Swarm::VelocityCounting() {
	double NewVelocity{ 0 }; //zmiana prêdkoœci cz¹stki
	for (int i = 0; i < Particles.size(); i++) {
		for (int j = 0; j < Particles[i].Velocity.size(); j++) {
			double r1 = Uniform_Mersenne(0, 1);
			double r2 = Uniform_Mersenne(0, 1);
			double r3 = Uniform_Mersenne(0, 1);

			NewVelocity += 0.2*r1*Particles[i].Velocity[j];  //czynnik - cz¹stka polega na sobie 
			NewVelocity += (0.2*r2*(Particles[i].ParticleBest[j] - Particles[i].Value)); //czynnik - cz¹stka polega na w³asnym doœwiadzceniu
			NewVelocity += (0.2*r3*(Particles[i].GlobalBest[j] - Particles[i].Value)); //czynnik - cz¹stka polega na doœwiadzceniu roju
			 Particles[i].Velocity[j] =  NewVelocity;
			 NewVelocity = 0;
		}
	}
}

void Swarm::PositionCounting() { //Przeliczenia pozycji cz¹stki
	for (int i = 0; i < Particles.size(); i++) {
		for (int j = 0; j < Particles[i].Position.size(); j++) {
			
			Particles[i].Position[j] += (Particles[SwarmSize-1].Position[j] - Particles[i].Position[j])*0.2;  //aktualizacja pozycji cz¹stki  
			
		}
	}
}

void Swarm::GlobalBestFinding() { //szukanie najlepszego rozwi¹zania dla ca³ego roju
	for (int i = 0; i < Particles.size(); i++) {
		for (int j = 0; j < Particles[i].Velocity.size(); j++) {
			if (Particles[i].GlobalBest[j] < Particles[SwarmSize-1].Value) {
				Particles[i].GlobalBest[j] = FunctionCalculating( Particles[SwarmSize - 1].Position[j]);    //przypisanie nowej globalnej najlepszej pozycji
			}
		}
	}
}

void Swarm::FindigPersonalBest() {   //szukanie najlepszego rzowi¹zania dla cz¹stki
	for (int i = 0; i < Particles.size(); i++) {
		for(int j = 0; j < Particles[i].Velocity.size(); j++) {
			if (Particles[i].Velocity[j] > Particles[i].ParticleBest[j]) {
				Particles[i].ParticleBest[j] = FunctionCalculating(Particles[i].Position[j]);   //przypisanie nowej najlepszej pozycji cz¹stki 
			}
		}
	}
}
void Swarm::RangeProtecting() {
	for (int i = 0; i < Particles.size(); i++) {
		for (int j = 0; j < Particles[i].Position.size(); j++) {
			if (Particles[i].Position[j] < LOW) {
				Particles[i].Position[j] = LOW;
			}
			else
				if (Particles[i].Position[j] > UP) {
					Particles[i].Position[j] = UP;
				}
		}
	}
}

void Swarm::SwarmCaluclating() {

	for (int k = 0; k < Particles.size(); k++) {
		for (int j = 0; j < Particles[k].Velocity.size(); j++) {

			Particles[k].Value = FunctionCalculating(Particles[k].Position);

		}
	}
	for (int p = 0; p < Particles.size(); p++) {
		for (int w = 0; w < Particles[p].Velocity.size(); w++) {

			Particles[p].Velocity[w] = Particles[SwarmSize - 1].Value - Particles[p].Value;

		}
	}
 }

void Swarm::ResultsPresenting(int iterator) {
	double Average{ 0 };
	for (int j = 0; j < Particles.size(); j++) {
		Average += Particles[j].Value;
	}
	
	cout << "| wartosc najlepszej czastki: " << Particles[SwarmSize - 1].Value << " | srednia wartosc roju: " << (Average / SwarmSize) << "| w iteracji: " << iterator+1 << "\n";

	Average = 0;

}

void Swarm::PSO() {//wywo³anie algorytmu

	SwarmGenerating();
	
		for (int i = 0; i < 50; i++) {
		
			SwarmCaluclating();

			Sorting();

            GlobalBestFinding();

			FindigPersonalBest();

			VelocityCounting();

			PositionCounting();

			RangeProtecting();

			ResultsPresenting(i);
				
		}
}