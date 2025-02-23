source("PfLOME_HealthState.R")
source("PfLOME_ImmuneState.R")

Human <- R6Class("Human",

                 ## Public Fields, Methods, and Initialization
                 public = list(

                   ## Initialization of Components

                   initialize = function(ixH = NA, age = NA, sex = NA, locH = NA, IncImm=T){
                     private$ixH = ixH
                     private$age = age
                     private$sex = sex
                     private$locH = locH
                     private$pathogen = Pathogen$new()
                     private$IncImm = IncImm
                     if(private$IncImm == T){
                        private$immuneState = ImmuneState$new()
                     }
                     private$healthState = HealthState$new()
                   },

                   ######## Infection Methods #########

                   infectHuman = function(t,pf){
                     if(private$IncImm == T){
                         BSImm = private$immuneState$get_BSImm()
                     }
                     else {
                       BSImm = 0
                     }
                     pfid = pf$get_pfid()
                     gtype = pf$get_gtype()
                     nptypes = pf$get_nptypes()
                     private$pathogen$add_Pf(t,pfid,gtype,BSImm,nptypes)
                     #typeImm = private$immuneState$get_typeImm(t,ptype)

                   },

                   ## write method to remove particular infection
                   ## **not currently used **
                   ## WAIT ARE PATHOGENS NOT CURRENTLY BEING CLEARED???!!!
                   clearPathogen = function(t, pfid){
                     private$pathogen$PfPathogen[[pfid]] = NULL
                     private$pathogen$set_PfMOI(private$pathogen$get_PfMOI()-1)
                   },

                   Treat = function(t,Drug){
                     private$healthState$Treat(t,Drug)
                   },

                   ########## Update Function #########

                   updateHuman = function(t,dt){
                     ##only update immunity if included
                     if(private$IncImm == T){
                        private$immuneState$update_immuneState(t,dt,self$get_Ptot())
                     }
                     private$healthState$update_healthState(t,dt,self$get_Ptot(),self$get_history()$RBC)
                     private$pathogen$update_pathogen(t,dt,self$get_PD())
                   },

                   ########## Accessors ##############

                   get_ixH = function(){
                     private$ixH
                   },

                   set_ixH = function(newixH){
                     private$ixH = newixH
                   },

                   get_age = function(){
                     private$age
                   },

                   set_age = function(newAge){
                     private$age = newAge
                   },

                   get_sex = function(){
                     private$sex
                   },

                   set_sex = function(newSex){
                     private$sex = newSex
                   },
                   get_immuneState = function(){
                     if(private$IncImm == T){
                        return(private$immuneState)
                     }
                     if(private$IncImm == F){
                       return(private$IncImm)
                     }
                   },

                   get_healthState = function(){
                     private$healthState
                   },

                   get_pathogen = function(){
                     private$pathogen
                   },

                   get_HRP2 = function(){
                     private$healthState$get_HRP2()
                   },

                   get_pLDH = function(){
                     private$healthState$get_pLDH()
                   },

                   get_Ptot = function(){
                     private$pathogen$get_Ptot()
                   },

                   get_Drug = function(){
                     private$healthState$get_Drug()
                   },

                   get_RxStart = function(){
                     private$healthState$get_RxStart()
                   },

                   get_PD = function(){
                     private$healthState$get_PD()
                   },

                   get_Fever = function(){
                     private$healthState$get_Fever()
                   },

                   get_PfMOI = function(){
                     private$pathogen$get_PfMOI()
                   },

                   get_history = function(){
                     if(private$IncImm==T){
                       c(private$pathogen$get_history(),
                       private$healthState$get_history(),
                       private$immuneState$get_history())
                     }
                     else {
                       c(private$pathogen$get_history(),
                       private$healthState$get_history())
                     }
                   },

                   get_IncImm = function(){
                     private$IncImm
                   }
                 ),

                 ## Private Fields
                 private = list(
                   ixH = NULL,
                   age = NULL,
                   sex = NULL,
                   locH = NULL,
                   IncPfPed = NULL,
                   pathogen = NULL,
                   IncImm = NULL,
                   immuneState = NULL,
                   healthState = NULL
                 )
)
