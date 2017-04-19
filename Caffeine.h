#ifndef _CAFFEINE_HEADER_
#define _CAFFEINE_HEADER_

// this is a base class of energy dependencies.
// it has sub-classes for each dependence

class caffeine {
	private:
		std::string dependence_name;
	public:
		caffeine(std::string label) : dependence_name(label) {};
       ~caffeine() {};

    	virtual std::string name() final { return dependence_name; };
    	virtual double sample(double E) = 0;  

};


class constant_dependence : public caffeine {
  private:
    double a;
  public:
     constant_dependence(std::string label, double A) :
     	caffeine(label), a(A) {};
    ~constant_dependence() {};

     double sample(double E) { return a; };
};

class inverse_sqrt_dependence : public caffeine {
  private:
    double a, b;
  public:
     inverse_sqrt_dependence(std::string label, double A, double B) :
     	caffeine(label), a(A), b(B) {};
    ~inverse_sqrt_dependence() {};

     double sample(double E) { return a + b / std::sqrt(E); };
};


#endif


