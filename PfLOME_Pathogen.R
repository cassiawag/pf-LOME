Pathogen <- R6Class("Pathogen",

                    public = list(

                      ## initialization of components

                      initialize = function(){
                        private$PfPathogen = list()
                        private$PfMOI = 0
                        private$Ptot = NaN
                        private$history = list()
                      },

                      ## add pf during infection

                      add_Pf = function(t,pfid,gtype,BSImm, nptypes){ # Where is type specific immunity?
                        pf = Pf$new(pfid, nptypes)
                        pf$set_gtype(gtype)
                        pf$set_PAR(pf$tentPAR(t,pfid))
                        pf$immuneMod_Tent(BSImm)
                        pf$set_Pt(pf$get_PAR()$MZ0)
                        pf$set_activeP(1)
                        ifelse(is.na(private$Ptot),self$set_Ptot(pf$get_PAR()$MZ0),self$set_Ptot(log10(10^pf$get_PAR()$MZ0+10^private$Ptot)))
                        private$PfPathogen[[pfid]] = pf
                        self$set_PfMOI(self$get_PfMOI()+1)
                      },

                      remove_Pf = function(t,pfid){
                        private$PfPathogen[[pfid]] = NULL ## null out object or just set active P/G to 0?
                        self$set_PfMOI(self$get_PfMOI()-1)
                      },

                      ## update pathogens

                      update_pathogen = function(t,dt,PD){
                        if (length(private$PfPathogen) > 0){
                          for(i in 1:length(private$PfPathogen)){
                            Pttemp = private$PfPathogen[[i]]$get_Pt()
                            private$PfPathogen[[i]]$update_Pf(t,dt,PD)
                            if(is.na(private$PfPathogen[[i]]$get_Pt()) & !is.na(Pttemp)){
                              self$set_PfMOI(private$PfMOI-1)
                            }
                          }
                          self$update_Ptot()
                        }
                        self$update_history()
                      },


                      ######### update methods ##########


                      update_Ptot = function(){
                        private$Ptot = NaN
                        for(i in 1:length(private$PfPathogen)){
                          private$Ptot = self$log10sum(c(private$Ptot, private$PfPathogen[[i]]$get_Pt()))
                        }
                      },

                      update_history = function(overwrite=F){
                        if(overwrite==F){
                          private$history$Ptot = c(private$history$Ptot,private$Ptot)
                          private$history$Gtot = c(private$history$Gtot,private$Gtot)
                          private$history$TE = c(private$history$TE,private$TE)
                          private$history$PfMOI = c(private$history$PfMOI,private$PfMOI)
                        }
                        if(overwrite==T){
                          private$history$Ptot[length(private$history$Ptot)] = private$Ptot
                          private$history$Gtot[length(private$history$Gtot)] = private$Gtot
                          private$history$TE[length(private$history$TE)] = private$TE
                          private$history$PfMOI[length(private$history$PfMOI)] = private$PfMOI
                        }
                      },

                      ######## accessors #########

                      get_Ptot = function(){
                        private$Ptot
                      },

                      set_Ptot = function(newPtot){
                        private$Ptot = newPtot
                      },

                      get_PfMOI = function(){
                        private$PfMOI
                      },

                      set_PfMOI = function(newPfMOI){
                        private$PfMOI = newPfMOI
                      },

                      get_Pf = function(){
                        private$PfPathogen
                      },

                      get_history = function(){
                        private$history
                      },

                      log10sum = function(x){
                        self$log10vals(log10(sum(10^self$log10vals(x), na.rm=TRUE)))
                      },

                      log10vals = function(x){
                        ifelse(!is.na(x) & is.finite(x) & x>=0, x, NaN)
                      }
                    ),


                    ########### private fields ##############


                    private = list(
                      PfPathogen = NULL,
                      Ptot = NULL,
                      Gtot = NULL,
                      PfMOI = NULL,
                      history = NULL
                    )

)

Pf <- R6Class("Pf",

              public = list(

                ## initialization of components

                initialize = function(pfid, nptypes){
                  private$pfid = pfid
                  private$Ptt = rep(NaN,10)
                  private$mnMaxPD = 10.5 # median peak parasitemia
                  private$mnPeakD = 20 # median peak time
                  private$mnMZ0 = 4.2 # median number of merozoites emerging from liver
                  private$mnDuration=200 # median infection duration
                  private$nptypes = nptypes
                  private$gtype = self$getGtype(nptypes)
                  private$ptype = self$getPtype(private$gtype,nptypes)
                },


                ######### setting g/p types


                getGtype = function(nptypes){
                  len <- length(nptypes)
                  gtype <- runif(len)
                  return(gtype)
                },

                getPtype = function(gtype,nptypes){
                  ptype = ceiling(gtype*nptypes)
                  return(ptype)
                },

                ############ accessors ############


                get_pfid = function(){
                  private$pfid
                },

                get_nptypes = function(){
                  private$nptypes
                },

                get_gtype = function(){
                  private$gtype
                },

                set_gtype = function(newgtype){
                  private$gtype = newgtype
                },
                get_ptype = function(){
                  private$ptype
                },

                set_ptype = function(newptype){
                  private$ptype = newptype
                },

                get_PAR = function(){
                  private$PAR
                },

                set_PAR = function(newPAR){
                  private$PAR = newPAR
                },

                get_Pt = function(){
                  private$Pt
                },

                set_Pt = function(newPt){
                  private$Pt = newPt
                },

                get_Ptt = function(){
                  private$Ptt
                },

                set_Ptt = function(newPtt){
                  private$Ptt = newPtt
                },

                get_activeP = function(){
                  private$activeP
                },

                set_activeP = function(newactiveP){
                  private$activeP = newactiveP
                },

                get_activeG = function(){
                  private$activeG
                },

                set_activeG = function(newactiveG){
                  private$activeG = newactiveG
                },

                ##average values for tent function, if not default
                set_mnMaxPD = function(mnMaxPD){
                  private$mnMaxPD = mnMaxPD
                  private$PAR = self$tentPAR(t=private$PAR$t0,pfid=private$pfid)
                },

                set_mnPeakD = function(mnPeakD){
                  private$mnPeakD = mnPeakD
                  private$PAR = self$tentPAR(t=private$PAR$t0,pfid=private$pfid)
                },

                set_mnMZ0 = function(mnMZ0){
                  private$mnMZ0 = mnMZ0
                  private$PAR = self$tentPAR(t=private$PAR$t0,pfid=private$pfid)
                },

                set_mnDuration = function(mnDuration){
                  private$mnDuration = mnDuration
                  private$PAR = self$tentPAR(t=private$PAR$t0,pfid=private$pfid)
                },

                ########## update methods ##########

                update_Pf = function(t,dt,PD){
                  for(i in 1:dt){
                    self$update_Pt(t,PD)
                    self$update_Ptt()
                  }
                },

                update_Pt = function(t,PD){
                  self$set_Pt(self$dPdt_tent(t,private$Pt,private$PAR,PD))
                  if(is.na(private$Pt)){
                    private$PAR$tEnd = t - private$PAR$t0
                    self$set_activeP(0)
                  }
                },

                update_Ptt = function(){
                  private$Ptt = self$shift(private$Ptt,1)
                  private$Ptt[1] = private$Pt
                },

                update_PD = function(PD){
                  private$PD = set_PD(PD)
                },


                ############### Tent Methods #################


                Pf.MaxPD = function(N=1, mn=private$mnMaxPD, vr=0.5){
                  #rnorm(N,mn,vr)
                  mn
                },

                Pf.PeakD = function(min=18,mn=private$mnPeakD){
                  #FIX STUB
                  # Day when parasitemia first peaks
                  #ceiling(min+rlnorm(1,log(3),.5))
                  min+(mn-min)
                },

                Pf.MZ0 = function(mn=private$mnMZ0){
                  #FIX STUB
                  #rnorm(1,4.2,.1)
                  mn
                },

                Pf.Duration = function(peakD,N=1,mn=private$mnDuration){
                  #FIX STUB
                  # Time to last parasitemia
                  #peakD + rgeom(N,1/mn)
                  peakD + mn
                },


                tentPAR = function(t,pfid){
                  mxPD          = self$Pf.MaxPD(mn=private$mnMaxPD)
                  peakD         = self$Pf.PeakD(mn=private$mnPeakD)
                  MZ0           = self$Pf.MZ0(mn=private$mnMZ0)
                  tEnd          = self$Pf.Duration(peakD,mn=private$mnDuration)

                  gr 		        = (mxPD-MZ0)/(peakD)
                  dr            = mxPD/(tEnd-peakD)
                  gtype         = private$gtype

                  list(
                    pfid	        = pfid,
                    t0  	        = t,
                    gr            = gr,
                    dr            = dr,
                    MZ0           = MZ0,
                    peakD         = peakD,
                    mxPD          = mxPD,
                    tEnd          = tEnd
                  )
                },

                gr_tent = function(t, PAR){with(PAR,{
                  ifelse(t<peakD, gr, -dr)
                })},

                dPdt_tent = function(t, P, PAR, PD=0, IM=0){with(PAR,{
                  age = ifelse(t>=t0, t-t0, 0)
                    P = ifelse(age>=1 & age<=tEnd,
                             pmin(mxPD, P + self$gr_tent(age,PAR))-PD-IM,
                             NaN)
                    ifelse(!is.na(P)&P>0, P, NaN)
                    return(P)
                })},

                shift = function(v,places,dir="right") {
                  places = places%%length(v)
                  if(places==0) return(v)
                  temp = rep(0,length(v))
                  if(dir=="left"){
                    places = length(v)-places
                    dir = "right"}
                  if(dir=="right"){
                    temp[1:places] = tail(v,places)
                    temp[(places+1):length(v)] = v[1:(length(v)-places)]
                    return(temp)
                  }
                },

                log10sum = function(x){
                  self$log10vals(log10(sum(10^self$log10vals(x), na.rm=TRUE)))
                },

                log10vals = function(x){
                  ifelse(!is.na(x) & is.finite(x) & x>=0, x, NaN)
                },

                immuneMod_Tent = function(BSImm){ ###### Don't have a way for type immunity to modulate this!!!
                  genImm              = 1-prod(1-BSImm) ##total strength of combined general immunity
                  ################### as a stub - have maximum effect of immunity with some effect inevitable (prop here 2/3 to 1/3)
                  p = 2/3
                  private$PAR$MZ0     = self$sigmoid01(genImm,.5,-1,private$PAR$MZ0*p)+private$PAR$MZ0*(1-p)
                  private$PAR$mxPD    = self$sigmoid01(genImm,.5,-1,private$PAR$mxPD*p)+private$PAR$mxPD*(1-p)
                  private$PAR$peakD   = self$sigmoid01(genImm,.5,-1,private$PAR$peakD*p)+private$PAR$peakD*(1-p)
                  private$PAR$tEnd    = self$sigmoid01(genImm,.5,-1,private$PAR$tEnd*p)+private$PAR$tEnd*(1-p)
                },

                sigmoid01 = function(x,xh,b,max) { ##defined on [0,1] --> [0,max]
                  #xh is 50th percentile in (0,1)
                  #b is slope param - b in R\{0}
                  a = tan(pi/2*xh)^b
                  max/(a/tan(pi/2*x)^b+1)
                }

              ),



              ############ private fields #################


              private = list(
                pfid = NULL,
                ## tentPars
                PAR = NULL,
                ## Parasite Densities (Pt = Asexual, Gt = Gametocyte,
                ## St = Sporozoite)
                activeP = NULL,
                activeG = NULL,
                activeS = NULL,
                Pt = NULL,
                Ptt = NULL,
                ## biological parameters
                gtype = NULL,
                ptype = NULL,
                nptypes = NULL,

                ##tent function dist'n pars
                mnMaxPD = NULL,
                mnPeakD = NULL,
                mnMZ0 = NULL,
                mnDuration = NULL
              )
)
