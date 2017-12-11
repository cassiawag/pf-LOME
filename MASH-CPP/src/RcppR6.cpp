// Generated by RcppR6 (0.2.4): do not edit by hand
#include <MASHcpp.h>

// [[Rcpp::export]]
MASHcpp::HumanEventQ HumanEventQ__ctor(int initQ) {
  return MASHcpp::HumanEventQ(initQ);
}
// [[Rcpp::export]]
Rcpp::List HumanEventQ__firstEvent(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_) {
  return obj_->firstEvent();
}
// [[Rcpp::export]]
double HumanEventQ__firstTime(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_) {
  return obj_->firstTime();
}
// [[Rcpp::export]]
void HumanEventQ__rmFirstEventFromQ(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_) {
  obj_->rmFirstEventFromQ();
}
// [[Rcpp::export]]
void HumanEventQ__rmTagFromQ(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_, std::string tag) {
  obj_->rmTagFromQ(tag);
}
// [[Rcpp::export]]
int HumanEventQ__get_queueN(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_) {
  return obj_->get_queueN();
}
// [[Rcpp::export]]
void HumanEventQ__addEvent2Q(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_, Rcpp::List event) {
  obj_->addEvent2Q(event);
}
// [[Rcpp::export]]
Rcpp::List HumanEventQ__get_EventQ(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_) {
  return obj_->get_EventQ();
}
// [[Rcpp::export]]
void HumanEventQ__clearQ(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ> obj_) {
  obj_->clearQ();
}

// [[Rcpp::export]]
MASHcpp::HistoryGeneric HistoryGeneric__ctor(int N) {
  return MASHcpp::HistoryGeneric(N);
}
// [[Rcpp::export]]
void HistoryGeneric__track_history(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryGeneric> obj_, double tEvent, std::string event) {
  obj_->track_history(tEvent, event);
}
// [[Rcpp::export]]
Rcpp::List HistoryGeneric__get_history(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryGeneric> obj_) {
  return obj_->get_history();
}

// [[Rcpp::export]]
MASHcpp::HistoryTravel HistoryTravel__ctor(int N) {
  return MASHcpp::HistoryTravel(N);
}
// [[Rcpp::export]]
void HistoryTravel__track_travel(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryTravel> obj_, double tTravel, int locationH) {
  obj_->track_travel(tTravel, locationH);
}
// [[Rcpp::export]]
Rcpp::List HistoryTravel__get_travelHistory(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryTravel> obj_) {
  return obj_->get_travelHistory();
}

// [[Rcpp::export]]
MASHcpp::humanPfSI humanPfSI__ctor(int PfID_init, double tInf_init, double b_init, double c_init, bool infected_init, bool chemoprophylaxis_init, int N) {
  return MASHcpp::humanPfSI(PfID_init, tInf_init, b_init, c_init, infected_init, chemoprophylaxis_init, N);
}
// [[Rcpp::export]]
std::vector<int> humanPfSI__get_PfID(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_PfID();
}
// [[Rcpp::export]]
void humanPfSI__push_PfID(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, int PfID_new) {
  obj_->push_PfID(PfID_new);
}
// [[Rcpp::export]]
int humanPfSI__back_PfID(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->back_PfID();
}
// [[Rcpp::export]]
std::vector<double> humanPfSI__get_tInf(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_tInf();
}
// [[Rcpp::export]]
void humanPfSI__push_tInf(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, double tInf_new) {
  obj_->push_tInf(tInf_new);
}
// [[Rcpp::export]]
double humanPfSI__get_b(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_b();
}
// [[Rcpp::export]]
void humanPfSI__set_b(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, double b_new) {
  obj_->set_b(b_new);
}
// [[Rcpp::export]]
double humanPfSI__get_c(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_c();
}
// [[Rcpp::export]]
void humanPfSI__set_c(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, double c_new) {
  obj_->set_c(c_new);
}
// [[Rcpp::export]]
std::vector<std::string> humanPfSI__get_vectorInf(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_vectorInf();
}
// [[Rcpp::export]]
void humanPfSI__push_vectorInf(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, std::string vectorInf_new) {
  obj_->push_vectorInf(vectorInf_new);
}
// [[Rcpp::export]]
bool humanPfSI__get_infected(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_infected();
}
// [[Rcpp::export]]
void humanPfSI__set_infected(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, bool infected_new) {
  obj_->set_infected(infected_new);
}
// [[Rcpp::export]]
bool humanPfSI__get_chemoprophylaxis(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_chemoprophylaxis();
}
// [[Rcpp::export]]
void humanPfSI__set_chemoprophylaxis(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, bool chemoprophylaxis_new) {
  obj_->set_chemoprophylaxis(chemoprophylaxis_new);
}
// [[Rcpp::export]]
void humanPfSI__track_history(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_, double tEvent, std::string event) {
  obj_->track_history(tEvent, event);
}
// [[Rcpp::export]]
Rcpp::List humanPfSI__get_history(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI> obj_) {
  return obj_->get_history();
}

// [[Rcpp::export]]
MASHcpp::mosquitoPfSI mosquitoPfSI__ctor(int PfID_init, std::string MosquitoID_init, double tInf_init, bool infected_init) {
  return MASHcpp::mosquitoPfSI(PfID_init, MosquitoID_init, tInf_init, infected_init);
}
// [[Rcpp::export]]
int mosquitoPfSI__get_PfID(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_) {
  return obj_->get_PfID();
}
// [[Rcpp::export]]
void mosquitoPfSI__set_PfID(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_, int PfID_new) {
  obj_->set_PfID(PfID_new);
}
// [[Rcpp::export]]
std::string mosquitoPfSI__get_MosquitoID(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_) {
  return obj_->get_MosquitoID();
}
// [[Rcpp::export]]
void mosquitoPfSI__set_MosquitoID(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_, std::string MosquitoID_new) {
  obj_->set_MosquitoID(MosquitoID_new);
}
// [[Rcpp::export]]
double mosquitoPfSI__get_tInf(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_) {
  return obj_->get_tInf();
}
// [[Rcpp::export]]
void mosquitoPfSI__set_tInf(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_, double tInf_new) {
  obj_->set_tInf(tInf_new);
}
// [[Rcpp::export]]
std::string mosquitoPfSI__get_humanInf(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_) {
  return obj_->get_humanInf();
}
// [[Rcpp::export]]
void mosquitoPfSI__set_humanInf(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_, std::string humanInf_new) {
  obj_->set_humanInf(humanInf_new);
}
// [[Rcpp::export]]
bool mosquitoPfSI__get_infected(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_) {
  return obj_->get_infected();
}
// [[Rcpp::export]]
void mosquitoPfSI__set_infected(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_, bool infected_new) {
  obj_->set_infected(infected_new);
}
// [[Rcpp::export]]
Rcpp::List mosquitoPfSI__get_history(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI> obj_) {
  return obj_->get_history();
}

// [[Rcpp::export]]
MASHcpp::humanPfMOI humanPfMOI__ctor(double b_init, double c_init, bool chemoprophylaxis_init) {
  return MASHcpp::humanPfMOI(b_init, c_init, chemoprophylaxis_init);
}
// [[Rcpp::export]]
std::vector<int> humanPfMOI__get_PfID(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  return obj_->get_PfID();
}
// [[Rcpp::export]]
std::vector<double> humanPfMOI__get_tInf(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  return obj_->get_tInf();
}
// [[Rcpp::export]]
int humanPfMOI__get_MOI(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  return obj_->get_MOI();
}
// [[Rcpp::export]]
double humanPfMOI__get_b(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  return obj_->get_b();
}
// [[Rcpp::export]]
void humanPfMOI__set_b(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_, double b_new) {
  obj_->set_b(b_new);
}
// [[Rcpp::export]]
double humanPfMOI__get_c(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  return obj_->get_c();
}
// [[Rcpp::export]]
void humanPfMOI__set_c(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_, double c_new) {
  obj_->set_c(c_new);
}
// [[Rcpp::export]]
bool humanPfMOI__get_chemoprophylaxis(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  return obj_->get_chemoprophylaxis();
}
// [[Rcpp::export]]
void humanPfMOI__set_chemoprophylaxis(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_, bool chemoprophylaxis_new) {
  obj_->set_chemoprophylaxis(chemoprophylaxis_new);
}
// [[Rcpp::export]]
void humanPfMOI__add_Infection(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_, int PfID_new, double tInf_new) {
  obj_->add_Infection(PfID_new, tInf_new);
}
// [[Rcpp::export]]
void humanPfMOI__clear_Infection(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_, int PfID_ix) {
  obj_->clear_Infection(PfID_ix);
}
// [[Rcpp::export]]
void humanPfMOI__clear_Infections(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  obj_->clear_Infections();
}
// [[Rcpp::export]]
std::vector<int> humanPfMOI__get_Infection(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI> obj_) {
  return obj_->get_Infection();
}

// [[Rcpp::export]]
MASHcpp::mosquitoPfMOI mosquitoPfMOI__ctor() {
  return MASHcpp::mosquitoPfMOI();
}
// [[Rcpp::export]]
std::vector<int> mosquitoPfMOI__get_PfID(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfMOI> obj_) {
  return obj_->get_PfID();
}
// [[Rcpp::export]]
int mosquitoPfMOI__get_MOI(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfMOI> obj_) {
  return obj_->get_MOI();
}
// [[Rcpp::export]]
void mosquitoPfMOI__add_infection(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfMOI> obj_, int PfID_new, double tInfected_new, double tInfectious_new) {
  obj_->add_infection(PfID_new, tInfected_new, tInfectious_new);
}
// [[Rcpp::export]]
std::vector<int> mosquitoPfMOI__get_infections(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfMOI> obj_, double tNow) {
  return obj_->get_infections(tNow);
}

// [[Rcpp::export]]
MASHcpp::RiskQ RiskQ__ctor() {
  return MASHcpp::RiskQ();
}
// [[Rcpp::export]]
int RiskQ__get_N(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  return obj_->get_N();
}
// [[Rcpp::export]]
void RiskQ__set_N(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, int N_new) {
  obj_->set_N(N_new);
}
// [[Rcpp::export]]
std::vector<std::string> RiskQ__get_who(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  return obj_->get_who();
}
// [[Rcpp::export]]
void RiskQ__push_who(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, std::string who_new) {
  obj_->push_who(who_new);
}
// [[Rcpp::export]]
std::vector<double> RiskQ__get_pTm(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  return obj_->get_pTm();
}
// [[Rcpp::export]]
void RiskQ__push_pTm(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, double pTm_new) {
  obj_->push_pTm(pTm_new);
}
// [[Rcpp::export]]
std::vector<double> RiskQ__get_w(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  return obj_->get_w();
}
// [[Rcpp::export]]
void RiskQ__push_w(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, double w_new) {
  obj_->push_w(w_new);
}
// [[Rcpp::export]]
void RiskQ__add_HumanHost(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, std::string who_new, double pTm_new, double w_new) {
  obj_->add_HumanHost(who_new, pTm_new, w_new);
}
// [[Rcpp::export]]
Rcpp::List RiskQ__get_HumanHost(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  return obj_->get_HumanHost();
}
// [[Rcpp::export]]
Rcpp::List RiskQ__get_HumanHostID(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, std::string ID) {
  return obj_->get_HumanHostID(ID);
}
// [[Rcpp::export]]
void RiskQ__clear_HumanHost(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  obj_->clear_HumanHost();
}
// [[Rcpp::export]]
int RiskQ__get_nOther(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  return obj_->get_nOther();
}
// [[Rcpp::export]]
void RiskQ__set_nOther(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, int nOther_new) {
  obj_->set_nOther(nOther_new);
}
// [[Rcpp::export]]
void RiskQ__add_OtherHost(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_, double otherW_new, std::string typeID_new) {
  obj_->add_OtherHost(otherW_new, typeID_new);
}
// [[Rcpp::export]]
Rcpp::List RiskQ__get_OtherHost(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ> obj_) {
  return obj_->get_OtherHost();
}

// [[Rcpp::export]]
MASHcpp::MatingQ MatingQ__ctor() {
  return MASHcpp::MatingQ();
}
// [[Rcpp::export]]
int MatingQ__get_N(MASHcpp::RcppR6::RcppR6<MASHcpp::MatingQ> obj_) {
  return obj_->get_N();
}
// [[Rcpp::export]]
void MatingQ__add_male2Q(MASHcpp::RcppR6::RcppR6<MASHcpp::MatingQ> obj_, std::string maleID_new, double mateFitness_new, int maleGenotype_new) {
  obj_->add_male2Q(maleID_new, mateFitness_new, maleGenotype_new);
}
// [[Rcpp::export]]
Rcpp::List MatingQ__get_MatingQ(MASHcpp::RcppR6::RcppR6<MASHcpp::MatingQ> obj_) {
  return obj_->get_MatingQ();
}
// [[Rcpp::export]]
void MatingQ__clear_MatingQ(MASHcpp::RcppR6::RcppR6<MASHcpp::MatingQ> obj_) {
  obj_->clear_MatingQ();
}

// [[Rcpp::export]]
MASHcpp::ImagoQ ImagoQ__ctor() {
  return MASHcpp::ImagoQ();
}
// [[Rcpp::export]]
void ImagoQ__clear_ImagoQ(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_) {
  obj_->clear_ImagoQ();
}
// [[Rcpp::export]]
void ImagoQ__clear_ImagoQTime(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_, double time) {
  obj_->clear_ImagoQTime(time);
}
// [[Rcpp::export]]
void ImagoQ__add_ImagoQ(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_, int N_new, double tEmerge_new, int genotype_new) {
  obj_->add_ImagoQ(N_new, tEmerge_new, genotype_new);
}
// [[Rcpp::export]]
double ImagoQ__track_ImagoQ(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_, double time) {
  return obj_->track_ImagoQ(time);
}
// [[Rcpp::export]]
int ImagoQ__get_N(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_) {
  return obj_->get_N();
}
// [[Rcpp::export]]
void ImagoQ__set_N(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_, int N_new) {
  obj_->set_N(N_new);
}
// [[Rcpp::export]]
Rcpp::List ImagoQ__get_ImagoQ(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_) {
  return obj_->get_ImagoQ();
}
// [[Rcpp::export]]
Rcpp::List ImagoQ__get_ImagoQTime(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ> obj_, double tNow, bool clear) {
  return obj_->get_ImagoQTime(tNow, clear);
}

// [[Rcpp::export]]
MASHcpp::EggQ EggQ__ctor() {
  return MASHcpp::EggQ();
}
// [[Rcpp::export]]
void EggQ__clear_EggQ(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_) {
  obj_->clear_EggQ();
}
// [[Rcpp::export]]
void EggQ__clear_EggQTime(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_, double time) {
  obj_->clear_EggQTime(time);
}
// [[Rcpp::export]]
void EggQ__add_EggQ(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_, int N_new, double tOviposit_new, int genotype_new) {
  obj_->add_EggQ(N_new, tOviposit_new, genotype_new);
}
// [[Rcpp::export]]
double EggQ__track_EggQ(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_, double time) {
  return obj_->track_EggQ(time);
}
// [[Rcpp::export]]
int EggQ__get_N(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_) {
  return obj_->get_N();
}
// [[Rcpp::export]]
void EggQ__set_N(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_, int N_new) {
  obj_->set_N(N_new);
}
// [[Rcpp::export]]
Rcpp::List EggQ__get_EggQ(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_) {
  return obj_->get_EggQ();
}
// [[Rcpp::export]]
Rcpp::List EggQ__get_EggQTime(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ> obj_, double tNow, bool clear) {
  return obj_->get_EggQTime(tNow, clear);
}

// [[Rcpp::export]]
MASHcpp::EL4P EL4P__ctor(int numGenotypes, double psi_new, double alpha_new, double p_new) {
  return MASHcpp::EL4P(numGenotypes, psi_new, alpha_new, p_new);
}
// [[Rcpp::export]]
void EL4P__oneStep(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  obj_->oneStep();
}
// [[Rcpp::export]]
void EL4P__oneStep_GEL4P(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double M, double eqAqua, double G, double lifespan) {
  obj_->oneStep_GEL4P(M, eqAqua, G, lifespan);
}
// [[Rcpp::export]]
void EL4P__burnIn_GEL4P(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double M, double eqAqua, double G, double lifespan, int tMax) {
  obj_->burnIn_GEL4P(M, eqAqua, G, lifespan, tMax);
}
// [[Rcpp::export]]
void EL4P__G2K_GEL4P(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double eqAqua, double G, double lifespan, int tMax) {
  obj_->G2K_GEL4P(eqAqua, G, lifespan, tMax);
}
// [[Rcpp::export]]
std::vector<double> EL4P__checkDX_GEL4P(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double eqAqua, double G, double lifespan, int tMax) {
  return obj_->checkDX_GEL4P(eqAqua, G, lifespan, tMax);
}
// [[Rcpp::export]]
void EL4P__addEggs(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double eggs_N, int genotype) {
  obj_->addEggs(eggs_N, genotype);
}
// [[Rcpp::export]]
Rcpp::List EL4P__get_allGenotypes(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  return obj_->get_allGenotypes();
}
// [[Rcpp::export]]
Rcpp::List EL4P__get_genotypeIx(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, int ix) {
  return obj_->get_genotypeIx(ix);
}
// [[Rcpp::export]]
double EL4P__get_psi(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  return obj_->get_psi();
}
// [[Rcpp::export]]
void EL4P__set_psi(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double psi_new) {
  obj_->set_psi(psi_new);
}
// [[Rcpp::export]]
double EL4P__get_alpha(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  return obj_->get_alpha();
}
// [[Rcpp::export]]
void EL4P__set_alpha(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double alpha_new) {
  obj_->set_alpha(alpha_new);
}
// [[Rcpp::export]]
double EL4P__get_p(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  return obj_->get_p();
}
// [[Rcpp::export]]
void EL4P__set_p(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, double p_new) {
  obj_->set_p(p_new);
}
// [[Rcpp::export]]
int EL4P__get_numGenotypes(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  return obj_->get_numGenotypes();
}
// [[Rcpp::export]]
double EL4P__get_totalLambda(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  return obj_->get_totalLambda();
}
// [[Rcpp::export]]
double EL4P__get_specificLambda(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, int ix) {
  return obj_->get_specificLambda(ix);
}
// [[Rcpp::export]]
void EL4P__reset(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_) {
  obj_->reset();
}
// [[Rcpp::export]]
void EL4P__set_pop(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P> obj_, Rcpp::List initPop) {
  obj_->set_pop(initPop);
}

// [[Rcpp::export]]
MASHcpp::ELP ELP__ctor(double alpha_new, double gamma_new, double psi_new, double tGrain_new) {
  return MASHcpp::ELP(alpha_new, gamma_new, psi_new, tGrain_new);
}
// [[Rcpp::export]]
void ELP__oneDay_aquaticDynamics(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double eggIn) {
  obj_->oneDay_aquaticDynamics(eggIn);
}
// [[Rcpp::export]]
double ELP__oneDay_Emergence(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  return obj_->oneDay_Emergence();
}
// [[Rcpp::export]]
std::vector<Rcpp::List> ELP__ecologicalSimulation(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double g, double f, double v, double L_init, double M_init, double tMax) {
  return obj_->ecologicalSimulation(g, f, v, L_init, M_init, tMax);
}
// [[Rcpp::export]]
Rcpp::List ELP__ecologicalSimulation2(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double g_init, double f_init, double v_init, double L_init, double M_init, double tMax) {
  return obj_->ecologicalSimulation2(g_init, f_init, v_init, L_init, M_init, tMax);
}
// [[Rcpp::export]]
double ELP__get_psi(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  return obj_->get_psi();
}
// [[Rcpp::export]]
void ELP__set_psi(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double psi_new) {
  obj_->set_psi(psi_new);
}
// [[Rcpp::export]]
double ELP__get_alpha(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  return obj_->get_alpha();
}
// [[Rcpp::export]]
void ELP__set_alpha(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double alpha_new) {
  obj_->set_alpha(alpha_new);
}
// [[Rcpp::export]]
double ELP__get_gamma(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  return obj_->get_gamma();
}
// [[Rcpp::export]]
void ELP__set_gamma(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double gamma_new) {
  obj_->set_gamma(gamma_new);
}
// [[Rcpp::export]]
double ELP__get_tGrain(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  return obj_->get_tGrain();
}
// [[Rcpp::export]]
Rcpp::List ELP__get_ELP(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  return obj_->get_ELP();
}
// [[Rcpp::export]]
void ELP__set_ELP(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double L1_new, double L2_new, double L3_new, double L4_new) {
  obj_->set_ELP(L1_new, L2_new, L3_new, L4_new);
}
// [[Rcpp::export]]
Rcpp::List ELP__get_parameters(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  return obj_->get_parameters();
}
// [[Rcpp::export]]
void ELP__set_parameters(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_, double alpha_new, double gamma_new, double psi_new, double tGrain_new) {
  obj_->set_parameters(alpha_new, gamma_new, psi_new, tGrain_new);
}
// [[Rcpp::export]]
void ELP__reset(MASHcpp::RcppR6::RcppR6<MASHcpp::ELP> obj_) {
  obj_->reset();
}

// [[Rcpp::export]]
MASHcpp::MosquitoFemaleHistory MosquitoFemaleHistory__ctor() {
  return MASHcpp::MosquitoFemaleHistory();
}
// [[Rcpp::export]]
void MosquitoFemaleHistory__historyInit(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory> obj_, Rcpp::Environment privateEnv) {
  obj_->historyInit(privateEnv);
}
// [[Rcpp::export]]
void MosquitoFemaleHistory__set_mateID(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory> obj_, std::string mateID_new) {
  obj_->set_mateID(mateID_new);
}
// [[Rcpp::export]]
void MosquitoFemaleHistory__historyTrack(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory> obj_, Rcpp::Environment privateEnv, bool alive) {
  obj_->historyTrack(privateEnv, alive);
}
// [[Rcpp::export]]
void MosquitoFemaleHistory__historyFeed(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory> obj_, Rcpp::Environment privateEnv) {
  obj_->historyFeed(privateEnv);
}
// [[Rcpp::export]]
void MosquitoFemaleHistory__calcBionomics(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory> obj_) {
  obj_->calcBionomics();
}
// [[Rcpp::export]]
Rcpp::List MosquitoFemaleHistory__exportHistory(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory> obj_) {
  return obj_->exportHistory();
}

// [[Rcpp::export]]
MASHcpp::MosquitoMaleHistory MosquitoMaleHistory__ctor() {
  return MASHcpp::MosquitoMaleHistory();
}
// [[Rcpp::export]]
void MosquitoMaleHistory__historyInit(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoMaleHistory> obj_, Rcpp::Environment privateEnv) {
  obj_->historyInit(privateEnv);
}
// [[Rcpp::export]]
void MosquitoMaleHistory__historyTrack(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoMaleHistory> obj_, Rcpp::Environment privateEnv, bool alive) {
  obj_->historyTrack(privateEnv, alive);
}
// [[Rcpp::export]]
Rcpp::List MosquitoMaleHistory__exportHistory(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoMaleHistory> obj_) {
  return obj_->exportHistory();
}


