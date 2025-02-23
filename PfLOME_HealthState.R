HealthState <- R6Class("HealthState",

                       public = list(

                         initialize = function(){
                           private$Fever = FALSE
                           private$feverThresh = 9
                           private$HRP2 = 0
                           private$RBC = 2.49
                           private$pLDH = 0
                           private$PD = 0
                           private$history = list()
                         },


                         ####### Accessors ########


                         get_Fever = function(){
                           private$Fever
                         },
                         get_feverThresh = function(){
                           private$feverThresh
                         },
                         set_feverThresh = function(feverThresh){
                           private$feverThresh = feverThresh
                         },
                         set_Fever = function(newFever){
                           private$Fever = newFever
                         },
                         get_HRP2 = function(){
                           private$HRP2
                         },
                         set_HRP2 = function(newHRP2){
                           private$HRP2 = newHRP2
                         },
                         get_pLDH = function(){
                           private$pLDH
                         },
                         set_pLDH = function(newpLDH){
                           private$pLDH = newpLDH
                         },
                         get_RBC = function(){
                           private$RBC
                         },
                         set_RBC = function(newRBC){
                           private$RBC = newRBC
                         },
                         get_history = function(){
                           private$history
                         },
                         get_RxStart = function(){
                           private$RxStart
                         },
                         get_Drug = function(){
                           private$Drug
                         },
                         get_PD = function(){
                           private$PD
                         },
                         get_RBChistory = function(){
                           private$history$RBC
                         },


                         ############ Update Methods ##############


                         update_healthState = function(t,dt,Ptot,RBCHist){
                           self$update_Fever(Ptot)
                           self$update_HRP2(Ptot)
                           self$update_pLDH(Ptot)
                           self$update_RBC(Ptot,RBCHist)
                           self$update_PD(t)
                           self$update_history()
                         },

                         update_Fever = function(Ptot){
                           if(!is.na(Ptot)){
                             private$Fever = ifelse(Ptot >= private$feverThresh, TRUE, FALSE)
                           }
                           if(is.na(Ptot)){
                             private$Fever = FALSE
                           }
                         },

                         update_HRP2 = function(Ptot){
                           a = .0019
                           b = log(2)/3.67
                           private$HRP2 = ifelse(is.na(Ptot),log10(10^private$HRP2-dt*b*10^private$HRP2),log10(10^private$HRP2+dt*(a*10^Ptot-b*10^private$HRP2)))
                         },

                         update_pLDH = function(Ptot){
                           a = .13
                           b = log(2)/2
                           c = private$pLDH
                           if(!is.na(Ptot)){
                            private$pLDH = log10(10^c + dt*(a*10^Ptot - b*10^c))
                           }else{
                            private$pLDH = log10((1-dt*b)*10^c)
                           }
                         },

                         update_RBC = function(Ptot,RBCHist){
                           a = log(2)/120 #RBC halflife
                           b = 1
                           c = 1.7
                           d = .5
                           e = 5*10^9
                           # rhat is now not a vector. If it should be, I definitely made an error with this bug fix.
                           if(t < 7){
                            rhat = 2.5
                           }else{
                            rhat = RBCHist[t-6]
                           }
                           r = private$RBC
                           if(is.nan(Ptot)){
                            private$RBC = r + dt*(- a*r + b*exp(-c*rhat))
                           }else{
                            private$RBC = r + dt*(- a*r + b*exp(-c*rhat) - d*10^Ptot/(e+10^Ptot)*r)
                           }
                         },

                         update_PD = function(t){
                           private$PD = self$getPD(t,private$RxStart,private$Drug)
                         },

                         update_history = function(){
                           private$history$Fever = c(private$history$Fever,private$Fever)
                           private$history$HRP2 = c(private$history$HRP2,private$HRP2)
                           private$history$pLDH = c(private$history$pLDH,private$pLDH)
                           private$history$RBC = c(private$history$RBC,private$RBC)
                           private$history$PD = c(private$history$PD,private$PD)
                         },


                         #################### Diagnostic Tests ####################


                         RDT = function(){
                           detect = 10
                           E1 = .1
                           E2 = .1
                           x = 10^private$HRP2
                           p = E1+(1-E1-E2)*self$sigmoidX(x,detect,3,13)
                           return(rbinom(1,1,p))
                         },

                         HSRDT = function(){
                           detect = 10
                           E1 = .1
                           E2 = .1
                           x = 10^private$HRP2
                           p = E1+(1-E1-E2)*self$sigmoidX(x,detect,3,13)
                           return(rbinom(1,1,p))
                         },

                         PCR = function(){
                           detect = 10
                           E1 = .1
                           E2 = .1
                           x = 10^private$HRP2
                           p = E1+(1-E1-E2)*self$sigmoidX(x,detect,3,13)
                           return(rbinom(1,1,p))
                         },

                         LAMP = function(){
                           detect = 10
                           E1 = .1
                           E2 = .1
                           x = 10^private$HRP2
                           p = E1+(1-E1-E2)*self$sigmoidX(x,detect,3,13)
                           return(rbinom(1,1,p))
                         },

                         LightMic = function(){
                           detect = 10
                           E1 = .1
                           E2 = .1
                           x = 10^private$HRP2
                           p = E1+(1-E1-E2)*self$sigmoidX(x,detect,3,13)
                           return(rbinom(1,1,p))
                         },

                         sigmoidX = function(X, X50=6, Xs=3, atMax=13){
                           pmin((1/(1+exp(-Xs*(X-X50))) - 1/(1+exp(Xs*X50)))/(1/(1+exp(-Xs*(atMax-X50))) - 1/(1+exp(Xs*X50))),1)
                         },


                         ####################### Rx methods #######################

                         Treat = function(t,Drug){
                           private$RxStart = c(private$RxStart,t)
                           private$Drug = c(private$Drug,1)
                         },

                         getPD = function(t, RxStart, Drug){
                           N = length(private$RxStart)
                           PD = 0
                           if(N > 0){
                            for(i in 1:N){
                              PDnew = self$PDi(t,private$RxStart[i],private$Drug[i])
                              if(PDnew>0){
                                PD = log10(10^PD+10^PDnew)
                              }
                            }
                           }
                           return(PD)
                         },

                         PDi = function(t, RxStart, Drug){
                           age = t-RxStart+1
                           PD = 0
                           if(age>=1 &  age<=RxRegister[[Drug]]$Duration){
                             PD = RxRegister[[Drug]]$PfPD[age]
                           }
                           return(PD)
                         }


                       ),


                       private = list(
                         Fever = NULL,
                         feverThresh = NULL,
                         HRP2 = NULL,
                         pLDH = NULL,
                         RBC = NULL,
                         RxStart = NULL,
                         Drug = NULL,
                         PD = NULL,
                         history = NULL
                       )

)
