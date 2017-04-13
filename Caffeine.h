#ifndef _CAFFEINE_HEADER_
#define _CAFFEINE_HEADER_

// this is a base class of energy dependencies.
// it has sub-classes for each dependence

class caffeine {
	private:
		
	protected:
		std::string dependence_name;
	public:
		caffeine() {};
       ~caffeine() {};

    	virtual std::string name() final { return dependence_name; };
    	virtual double sample(double E) = 0;  

};


class constant_dependence : public caffeine {
  private:
    
  public:
     constant_dependence()
     	{ dependence_name = "constant"; };
    ~constant_dependence() {};

     double sample(double E);
};

class inverse_sqrt_dependence : public caffeine {
  private:
    
  public:
     inverse_sqrt_dependence()
     	{ dependence_name = "inverse sqrt"; };
    ~inverse_sqrt_dependence() {};

     double sample(double E);
};


#endif

