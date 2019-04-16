/* ################################################################################
 *      __  ______   __________  ____
 *     /  |/  /   | / ____/ __ \/ __ \
 *    / /|_/ / /| |/ /   / /_/ / / / /
 *   / /  / / ___ / /___/ _, _/ /_/ /
 *  /_/  /_/_/  |_\____/_/ |_|\____/
 *
 *  A very simple version for the PRISM data
 *  Stand-alone PfSI with no mosquitos (EIR is forcing)
 *
################################################################################ */

/* Rcpp bits */
#include <Rcpp.h>

// [[Rcpp::plugins(cpp14)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

/* STL includes */
#include <string>
#include <functional>
#include <algorithm>
#include <memory>


/* global simulation time */
static unsigned int tnow_global = 0;

/* global parameters */
static double DurationPf = 200.0;
static double b = 0.55;


/* ################################################################################
* sample waiting time (hazard functions for events)
################################################################################ */

/* duration of infection */
double pfsi_ttClearPf(){
 double recover = R::rexp(DurationPf);
 // std::cout << "ttClearPf random number: " << recover << "\n";
 return recover;
};

/* no latent period */
double psfi_ttInfectionPf(){
  return 0;
}

/* ################################################################################
 * generic event class (abstract base)
################################################################################ */

class event {
public:

  /* constructor */
  event(std::string tag_, double tEvent_, std::function<void()> eventF_):
    tag(tag_),tEvent(tEvent_),eventF(eventF_) {};

  /* destructor */
  virtual ~event(){};

  /* move operators */
  event(event&&) = default;
  event& operator=(event&&) = default;

  /* copy operators */
  event(event&) = default;
  event& operator=(event&) = default;

  /* print (debugging) */
  void print(){
    std::cout << "event -- tag: " << tag << ", tEvent: " << tEvent << std::endl;
  };

  /* comparison for sorting */
  bool operator<(event e) const {
    return tEvent < e.tEvent;
  };

  /* information for event */
  std::string                        tag;
  double                             tEvent;
  std::function<void()>              eventF;

};


/* ################################################################################
* PfSI event declarations
################################################################################ */

/* need to forward declare events here */
class human;

/* simulated infectious bite; tag: PfSI_SimBite */
class e_pfsi_bite : public event {
public:
  /* constructor */
  e_pfsi_bite(double tEvent_, human* h);

  /* destructor */
  ~e_pfsi_bite(){};

};

/* start a PfSI infection; tag: PfSI_infection */
class e_pfsi_infect : public event {
public:
  /* constructor */
  e_pfsi_infect(double tEvent_, human* h);

  /* destructor */
  ~e_pfsi_infect(){};
};


/* end a PfSI infection; tag: PfSI_recovery */
class e_pfsi_recover : public event {
public:
  /* constructor */
  e_pfsi_recover(double tEvent_, human* h);

  /* destructor */
  ~e_pfsi_recover(){};
};


/* ################################################################################
* human class
################################################################################ */

using eventP = std::unique_ptr<event>;

/* the class definition */
class human {
public:

  /* constructor & destructor */
  human(const int id_, const std::vector<double>& EIR_size_, const std::vector<double>& EIR_prob_, const std::string init) :
    id(id_), tnow(0.0), state("S"), EIR_size(EIR_size_), EIR_prob(EIR_prob_), bites(0) {
      if(init.compare("I") == 0){
        addEvent2Q(e_pfsi_infect(0.0,this));
      }
    };
  ~human(){};

  /* move operators */
  human(human&&) = default;
  human& operator=(human&&) = default;

  /* copy operators */
  human(human&) = delete;
  human& operator=(human&) = delete;

  /* print */
  void print();

  /* accessors */
  u_int                 get_id(){return id;};
  double                get_tnow(){return tnow;};

  std::string           get_state(){return state;};
  void                  set_state(const std::string s){state = s;};
  int                   get_bites(){return bites;};


  /* event queue related functions */
  void                  addEvent2Q(event&& e);
  void                  rmTagFromQ(const std::string &tag);
  void                  fireEvent();

  /* interface */
  void                  simulate();

private:

  /* basic fields */
  u_int                 id; /* my id */
  double                tnow; /* my local simulation time (time of last jump) */
  std::string           state;

  std::vector<eventP>   eventQ;

  /* biting: for nbinom(size,prob) parameterization */
  std::vector<double>   EIR_size;
  std::vector<double>   EIR_prob;

  /* history */
  int                   bites;

  /* called by simulate */
  void                  queue_bites();

};




/* ################################################################################
 * event queue
################################################################################ */

/* add an event to my queue */
void human::addEvent2Q(event&& e){
  eventQ.emplace_back(std::make_unique<event>(e));
  std::sort(eventQ.begin(),eventQ.end(),[](const std::unique_ptr<event>& e1, const std::unique_ptr<event>& e2){
    return e1->tEvent < e2->tEvent;
  });
};

/* remove future queued events */
void human::rmTagFromQ(const std::string &tag){
  eventQ.erase(std::remove_if(eventQ.begin(),eventQ.end(),
                              [tag](eventP& e){
                                return(tag.compare(e->tag)==0);
                              }),eventQ.end());
};

/* fire the first event */
void human::fireEvent(){
  if(eventQ.size() > 0){
    tnow = eventQ.front()->tEvent; /* update local simulation time */
    eventQ.front()->eventF();
    eventQ.erase(eventQ.begin());
  }
};


/* ################################################################################
 * PfSI events
################################################################################ */

/* simulated biting event */
e_pfsi_bite::e_pfsi_bite(double tEvent_, human* h):
  event("PfSI_SimBite",tEvent_,[tEvent_,h](){

    // std::cout << "simbite occuring to " << h->get_id() << "\n";

    /* transmission efficiency */
    if(R::runif(0.0, 1.0) < b){
      double tInfStart = tEvent_ + psfi_ttInfectionPf();
      h->addEvent2Q(e_pfsi_infect(tInfStart,h));
    }

  })
{};

/* infection */
e_pfsi_infect::e_pfsi_infect(double tEvent_, human* h):
  event("PfSI_infection",tEvent_,[tEvent_,h](){

    /* no superinfection, and chx blocks new infections */
    if(h->get_state().compare("S") == 0){

      // std::cout << "infection occuring to " << h->get_id() << " at time " << tnow_global << "\n";

      /* i'm now infected */
      h->set_state("I");

      // std::cout << "checking my state  " << h->get_state() << "\n";

      /* queue clearance event */
      double tEnd = tEvent_ + pfsi_ttClearPf();

      // std::cout << "it will clear at  " <<tEnd<< "\n";

      h->addEvent2Q(e_pfsi_recover(tEnd,h));

    }

  })
{};

/* recovery */
e_pfsi_recover::e_pfsi_recover(double tEvent_, human* h):
 event("PfSI_recovery",tEvent_,[tEvent_,h](){

   /* i can only recover if i'm infected */
   if(h->get_state().compare("I") == 0){
     h->set_state("S");
   }

 })
{};


/* ################################################################################
 * PfSI simulation loop
################################################################################ */

void human::simulate(){

  // std::cout << "simulating:  " << id << "\n";

  /* fire all events that occur on this time step */
  while(eventQ.size() > 0 && eventQ.front()->tEvent < tnow_global){
    fireEvent();
  }

  /* queue bites */
  queue_bites();

};


void human::queue_bites(){

  // std::cout << "queue_bites occuring to " << id << "\n";

  /* parameters of nbinom biting */
  double size = EIR_size.at(tnow_global);
  double prob = EIR_prob.at(tnow_global);

  // std::cout << "size: " << size << " prob: " << prob << "\n";


  bites = (int)R::rnbinom(size, prob);

  // std::cout << "they got  " << bites << " many bites\n";

  if(bites > 0){
    for(size_t i=0; i<bites; i++){
      addEvent2Q(e_pfsi_bite(tnow_global,this));
    }
  }

};


/* ################################################################################
 * Run a simulation from R
################################################################################ */

/* human pointer */
using humanP = std::unique_ptr<human>;

// [[Rcpp::export]]
Rcpp::List tiny_pfsi(const unsigned int tmax,
                     const size_t nh,
                     const Rcpp::StringVector init,
                     const Rcpp::List EIR_size,
                     const Rcpp::List EIR_prob,
                     const bool pb){

  /* checks */
  if(nh != EIR_size.size() || nh != EIR_prob.size()){
    Rcpp::stop("number of humans to simulate must be same as length of EIR size and prob lists");
  }

  tnow_global = 0;

  /* set up our ensemble of people */
  std::vector<humanP> humans;
  humans.reserve(nh);

  for(size_t i=0; i<nh; i++){

    humans.emplace_back(std::make_unique<human>(
      i,
      Rcpp::as<std::vector<double> >(EIR_size.at(i)),
      Rcpp::as<std::vector<double> >(EIR_prob.at(i)),
      Rcpp::as< std::string >(init(i))
    ));

  }

  /* output matrix */
  Rcpp::IntegerMatrix out_mat(tmax,2);
  Rcpp::IntegerMatrix out_bites(tmax,nh);

  /* track progress */
  Progress progbar(tmax, pb);

  /* run simulation */
  while(tnow_global < tmax){

    if(Progress::check_abort()){
      Rcpp::stop("user abort detected");
    };

    /* sim humans */
    for(auto& h : humans){
      h->simulate();
    }

    /* write output */
    for(auto& h : humans){

      /* write bites */
      out_bites.at(tnow_global,h->get_id()) = h->get_bites();

      /* write states */
      // std::cout << "writing states .. " << h->get_state() << "\n";
      // std::string mystate = h->get_state();
      // std::cout << "mystate: " << mystate  << " i'm " << h->get_id() << " and its time: " << tnow_global << "\n";
      if(h->get_state().compare("S") == 0){
        out_mat.at(tnow_global,0) += 1;
      } else {
        // std::cout << "INFECTED\n";
        out_mat.at(tnow_global,1) += 1;
      }
    }

    progbar.increment();
    tnow_global++;
  }

  // // debug
  // std::vector<double> test = Rcpp::as<std::vector<double> >(EIR_size.at(0));
  // for(auto it : test){
  //   std::cout << it << " - ";
  // }

  /* return a list */
  return Rcpp::List::create(Rcpp::Named("states") = out_mat,
                            Rcpp::Named("bites") = out_bites);
};
