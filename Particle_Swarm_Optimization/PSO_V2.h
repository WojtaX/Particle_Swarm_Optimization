class Particle {
public:
	vector<double>Position; 
	vector<double> Velocity;
	vector<double> ParticleBest;
	vector<double> GlobalBest{ 0 };
	double LOW{ -100 };
	double UP{ 100 };
	double Value{ 0 };
	int Dimension{ 1 };

	Particle() : Position{ 0 }, Velocity{ 0 }, ParticleBest{ 0 } {}
	Particle(int Dimension);
};

class Swarm {
public:
	vector<Particle> Particles;
	
	int SwarmSize{ 5 };
	double LOW{ 0 };
	double UP{ 20 };

	int Dimension{ 1 };
	void Sorting();//vector<Particle>
	double FunctionCalculating(vector<double>Position);  //funkjca do optymalizacji
	double FunctionCalculating(double Position);
	void SwarmGenerating();
	void VelocityCounting(); 
	void PositionCounting();
	void GlobalBestFinding();
	void FindigPersonalBest();
	void RangeProtecting();
	void SwarmCaluclating();
	void ResultsPresenting(int iterator);
	void PSO(); //Funkcja wywo³uj¹ca wszytkie pozosta³e 
 };