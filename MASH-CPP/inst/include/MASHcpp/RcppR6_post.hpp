// Generated by RcppR6 (0.2.4): do not edit by hand
#ifndef _MASHCPP_RCPPR6_POST_HPP_
#define _MASHCPP_RCPPR6_POST_HPP_

#include <Rcpp.h>
#include <MASHcpp/RcppR6_support.hpp>

namespace MASHcpp {
namespace RcppR6 {
namespace traits {
template <> inline std::string   class_name_r<MASHcpp::HumanEventQ >() {return "HumanEventQ";}
template <> inline std::string   package_name<MASHcpp::HumanEventQ >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::HumanEventQ >() {return ".R6_HumanEventQ";}
template <> inline std::string   class_name_r<MASHcpp::HistoryGeneric >() {return "HistoryGeneric";}
template <> inline std::string   package_name<MASHcpp::HistoryGeneric >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::HistoryGeneric >() {return ".R6_HistoryGeneric";}
template <> inline std::string   class_name_r<MASHcpp::HistoryTravel >() {return "HistoryTravel";}
template <> inline std::string   package_name<MASHcpp::HistoryTravel >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::HistoryTravel >() {return ".R6_HistoryTravel";}
template <> inline std::string   class_name_r<MASHcpp::humanPfSI >() {return "humanPfSI";}
template <> inline std::string   package_name<MASHcpp::humanPfSI >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::humanPfSI >() {return ".R6_humanPfSI";}
template <> inline std::string   class_name_r<MASHcpp::mosquitoPfSI >() {return "mosquitoPfSI";}
template <> inline std::string   package_name<MASHcpp::mosquitoPfSI >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::mosquitoPfSI >() {return ".R6_mosquitoPfSI";}
template <> inline std::string   class_name_r<MASHcpp::humanPfMOI >() {return "humanPfMOI";}
template <> inline std::string   package_name<MASHcpp::humanPfMOI >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::humanPfMOI >() {return ".R6_humanPfMOI";}
template <> inline std::string   class_name_r<MASHcpp::mosquitoPfMOI >() {return "mosquitoPfMOI";}
template <> inline std::string   package_name<MASHcpp::mosquitoPfMOI >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::mosquitoPfMOI >() {return ".R6_mosquitoPfMOI";}
template <> inline std::string   class_name_r<MASHcpp::RiskQ >() {return "RiskQ";}
template <> inline std::string   package_name<MASHcpp::RiskQ >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::RiskQ >() {return ".R6_RiskQ";}
template <> inline std::string   class_name_r<MASHcpp::ImagoQ >() {return "ImagoQ";}
template <> inline std::string   package_name<MASHcpp::ImagoQ >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::ImagoQ >() {return ".R6_ImagoQ";}
template <> inline std::string   class_name_r<MASHcpp::EggQ >() {return "EggQ";}
template <> inline std::string   package_name<MASHcpp::EggQ >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::EggQ >() {return ".R6_EggQ";}
template <> inline std::string   class_name_r<MASHcpp::EL4P >() {return "EL4P";}
template <> inline std::string   package_name<MASHcpp::EL4P >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::EL4P >() {return ".R6_EL4P";}
template <> inline std::string   class_name_r<MASHcpp::ELPool >() {return "ELPool";}
template <> inline std::string   package_name<MASHcpp::ELPool >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::ELPool >() {return ".R6_ELPool";}
template <> inline std::string   class_name_r<MASHcpp::MosquitoFemaleHistory >() {return "MosquitoFemaleHistory";}
template <> inline std::string   package_name<MASHcpp::MosquitoFemaleHistory >() {return "MASHcpp";}
template <> inline std::string generator_name<MASHcpp::MosquitoFemaleHistory >() {return ".R6_MosquitoFemaleHistory";}
}
}
}

namespace Rcpp {
template <typename T>
SEXP wrap(const MASHcpp::RcppR6::RcppR6<T>& x) {
  return x.to_R6();
}

namespace traits {
template <typename T>
class Exporter<MASHcpp::RcppR6::RcppR6<T> > {
public:
  Exporter(SEXP x) : obj(MASHcpp::RcppR6::RcppR6<T>(x)) {}
  inline MASHcpp::RcppR6::RcppR6<T> get() { return obj; }
private:
  MASHcpp::RcppR6::RcppR6<T> obj;
};
}

template <> inline SEXP wrap(const MASHcpp::HumanEventQ& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ>(x));
}
template <> inline MASHcpp::HumanEventQ as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::HumanEventQ>(x));
}
template <> inline SEXP wrap(const MASHcpp::HistoryGeneric& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryGeneric>(x));
}
template <> inline MASHcpp::HistoryGeneric as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryGeneric>(x));
}
template <> inline SEXP wrap(const MASHcpp::HistoryTravel& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryTravel>(x));
}
template <> inline MASHcpp::HistoryTravel as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::HistoryTravel>(x));
}
template <> inline SEXP wrap(const MASHcpp::humanPfSI& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI>(x));
}
template <> inline MASHcpp::humanPfSI as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfSI>(x));
}
template <> inline SEXP wrap(const MASHcpp::mosquitoPfSI& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI>(x));
}
template <> inline MASHcpp::mosquitoPfSI as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfSI>(x));
}
template <> inline SEXP wrap(const MASHcpp::humanPfMOI& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI>(x));
}
template <> inline MASHcpp::humanPfMOI as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::humanPfMOI>(x));
}
template <> inline SEXP wrap(const MASHcpp::mosquitoPfMOI& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfMOI>(x));
}
template <> inline MASHcpp::mosquitoPfMOI as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::mosquitoPfMOI>(x));
}
template <> inline SEXP wrap(const MASHcpp::RiskQ& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ>(x));
}
template <> inline MASHcpp::RiskQ as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::RiskQ>(x));
}
template <> inline SEXP wrap(const MASHcpp::ImagoQ& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ>(x));
}
template <> inline MASHcpp::ImagoQ as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::ImagoQ>(x));
}
template <> inline SEXP wrap(const MASHcpp::EggQ& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ>(x));
}
template <> inline MASHcpp::EggQ as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::EggQ>(x));
}
template <> inline SEXP wrap(const MASHcpp::EL4P& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P>(x));
}
template <> inline MASHcpp::EL4P as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::EL4P>(x));
}
template <> inline SEXP wrap(const MASHcpp::ELPool& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::ELPool>(x));
}
template <> inline MASHcpp::ELPool as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::ELPool>(x));
}
template <> inline SEXP wrap(const MASHcpp::MosquitoFemaleHistory& x) {
  return wrap(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory>(x));
}
template <> inline MASHcpp::MosquitoFemaleHistory as(SEXP x) {
  return *(MASHcpp::RcppR6::RcppR6<MASHcpp::MosquitoFemaleHistory>(x));
}
}

#endif
