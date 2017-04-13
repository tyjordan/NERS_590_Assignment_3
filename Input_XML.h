#ifndef _INPUT_XML_HEADER_
#define _INPUT_XML_HEADER_


void Input_Problem_Data
( unsigned long long *Nsamples , bool *continuous_eng , bool *time_tracking , bool *split_roulette , 
  std::vector< std::shared_ptr< distribution<double> > > *double_distributions ,
  std::vector< std::shared_ptr< distribution<int>    > > *int_distributions ,
  std::vector< std::shared_ptr< distribution<point>  > > *point_distributions ,
  std::vector< std::shared_ptr< caffeine > > *eng_dependences ,
  std::vector< std::shared_ptr<nuclide> > *nuclides ,
  std::vector< std::shared_ptr<material> > *materials ,
  std::vector< std::shared_ptr< surface > > *surfaces ,
  std::vector< std::shared_ptr< cell > > *cells ,
  std::vector< std::shared_ptr< estimator > > *estimators ,
  std::shared_ptr< source > *src);

#endif
